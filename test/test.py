import os
import pytest
import docker
import csv
import json
import hashlib
import tarfile

from datetime import datetime
from reference_directory import ReferenceDirectory

dirname = os.path.dirname(__file__)
with open(os.path.join(dirname, 'test_data.json')) as config:
    TEST_DATA_CONFIG = json.load(config)

HOME_DIRECTORY = os.path.expanduser('~')

GRAPE_DOCKERFILE = 'containers/snakemake/Dockerfile'
GRAPE_IMAGE_TAG = 'genx_relatives:latest'

REFERENCE_DIRECTORY = os.path.join(HOME_DIRECTORY, 'ref_test')
CONTAINER_REFERENCE_DIRECTORY = '/media/ref'
CONTAINER_WORKING_DIRECTORY = '/media/data'
TESTING_REAL_DATA_DIRECTORY = os.path.join(HOME_DIRECTORY, 'test_data')
CONTAINER_TESTING_REAL_DATA_DIRECTORY = '/media/test_data'
KHAZAR_VCF = os.path.join(CONTAINER_TESTING_REAL_DATA_DIRECTORY, 'khazar314.vcf.gz')
AADR_VCF = os.path.join(CONTAINER_TESTING_REAL_DATA_DIRECTORY, 'aadr.reheaded.vcf.gz')
METRICS_FILEPATH = 'results/metrics.tsv'
RELATIVES_FILEPATH = 'results/relatives.tsv'
AADR_SAMPLES_FILEPATH = os.path.join(TESTING_REAL_DATA_DIRECTORY, 'aadr_samples.csv')


def _download_test_data(test_data_directory):
    url = TEST_DATA_CONFIG['download']['url']
    key = TEST_DATA_CONFIG['download']['azure_public_key']
    file = TEST_DATA_CONFIG['download']['file']
    md5_sum = TEST_DATA_CONFIG['download']['md5']
    tar = os.path.join(test_data_directory, file)

    os.system(f'wget "{url}{key}" -O {tar} --tries 50')

    with open(tar, 'rb') as test_data_tar:
        data = test_data_tar.read()
        md5_returned = hashlib.md5(data).hexdigest()
    if md5_returned == md5_sum:
        test_tar = tarfile.open(tar)
        test_tar.extractall(test_data_directory)
        test_tar.close()
        os.remove(tar)
        return True
    else:
        return False


def _read_metrics_file(filepath):
    metrics = {}
    with open(filepath, 'r') as metrics_file:
        reader = csv.DictReader(metrics_file, delimiter='\t')
        for row in reader:
            degree = row['True Degree']
            metrics[degree] = {
                'Precision': float(row['Precision']),
                'Recall': float(row['Recall'])
            }

    return metrics


def _read_samples_file(filepath):
    samples = {}
    with open(filepath, 'r') as samples_file:
        reader = csv.DictReader(samples_file)
        for row in reader:
            sample = row['id']
            samples[sample] = row['date']

    return samples


def _read_relatives_file(filepath):
    relatives = []
    with open(filepath, 'r') as relatives_file:
        reader = csv.DictReader(relatives_file, delimiter='\t')
        for row in reader:
            relative = (row['id1'], row['id2'])
            relatives.append(relative)

    return relatives


@pytest.fixture
def docker_client():
    client = docker.from_env()
    return client


@pytest.fixture
def grape_image(docker_client):
    """
    Build Docker image to evaluate tests.
    """
    docker_client.images.build(
        path='.', dockerfile=GRAPE_DOCKERFILE, tag=GRAPE_IMAGE_TAG,
        rm=True, container_limits={'memory': 8 * 1024 * 1024 * 1024}
    )
    print(f'\nCreated docker image with tag {GRAPE_IMAGE_TAG}')
    yield docker_client.images.get(GRAPE_IMAGE_TAG)

    # Fixture teardown to remove GRAPE Docker image
    docker_client.images.remove(GRAPE_IMAGE_TAG, force=True, noprune=False)
    print(f'\nRemoved docker image with tag {GRAPE_IMAGE_TAG}')


@pytest.fixture
def test_data_directory():
    content = TEST_DATA_CONFIG['content']
    actual_content = {}

    if not os.path.exists(TESTING_REAL_DATA_DIRECTORY):
        os.makedirs(TESTING_REAL_DATA_DIRECTORY)
        print('\nNo test data found, new test data archive will be downloaded!')
        if not _download_test_data(TESTING_REAL_DATA_DIRECTORY):
            raise Exception('Test data archive download failed!')

    for root, _, filenames in os.walk(TESTING_REAL_DATA_DIRECTORY):
        for filename in filenames:
            filepath = os.path.join(root, filename)
            relative_path = os.path.relpath(filepath, TESTING_REAL_DATA_DIRECTORY)
            actual_content[relative_path] = os.path.getsize(filepath)

    if actual_content != content:
        print('\nCurrent test data files does not match "test_data.json", new test data archive will be downloaded!')
        if not _download_test_data(TESTING_REAL_DATA_DIRECTORY):
            raise Exception('Test data archive download failed!')

    return TESTING_REAL_DATA_DIRECTORY


@pytest.fixture
def reference_directory(docker_client, grape_image) -> ReferenceDirectory:
    reference_directory = ReferenceDirectory(REFERENCE_DIRECTORY)

    if not reference_directory.is_valid():
        print('\nCurrrent reference data seem not to match '
              '"reference_directory_content.json", new reference data archive will be downloaded!')
        command = f'launcher.py reference --use-bundle --ref-directory {CONTAINER_REFERENCE_DIRECTORY} ' \
            '--phase --impute --real-run'
        volumes = {
            reference_directory.path: {'bind': CONTAINER_REFERENCE_DIRECTORY, 'mode': 'rw'}
        }

        docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=command, volumes=volumes)

    return reference_directory


@pytest.fixture(scope='function')
def working_directory(request):
    utc_timestamp = datetime.utcnow().strftime("%Y%m%d-%H%M%S-utc")
    working_directory_name = '-'.join([request.param, utc_timestamp])
    working_directory_path = os.path.join(HOME_DIRECTORY, working_directory_name)

    yield working_directory_path

    # Fixture teardown to remove working directory
    shutil.rmtree(working_directory_path)


@pytest.mark.parametrize('working_directory', ['ibis'], indirect=True)
def test_simulation_ibis(docker_client, grape_image, reference_directory, working_directory):
    volumes = {
        reference_directory.path: {'bind': CONTAINER_REFERENCE_DIRECTORY, 'mode': 'ro'},
        working_directory: {'bind': CONTAINER_WORKING_DIRECTORY, 'mode': 'rw'}
    }

    simulate_command = f'launcher.py simulate --ref-directory {CONTAINER_REFERENCE_DIRECTORY} --cores 8 '\
                       f'--directory {CONTAINER_WORKING_DIRECTORY} --flow ibis --assembly hg37 --seed 42 --real-run' \

    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=simulate_command, volumes=volumes)

    # Read file result with simulation metrics
    metrics = _read_metrics_file(os.path.join(working_directory, METRICS_FILEPATH))

    # Validate simultation metrics
    assert metrics['1']['Recall'] > 0.99 and metrics['1']['Precision'] > 0.99
    assert metrics['2']['Recall'] > 0.99 and metrics['2']['Precision'] > 0.99
    assert metrics['3']['Recall'] > 0.99 and metrics['3']['Precision'] > 0.99

    assert metrics['4']['Recall'] > 0.90 and metrics['4']['Precision'] > 0.95
    assert metrics['5']['Recall'] > 0.90 and metrics['5']['Precision'] > 0.95
    assert metrics['6']['Recall'] > 0.80 and metrics['6']['Precision'] > 0.90

    assert metrics['7']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['8']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['9']['Recall'] > 0 and metrics['1']['Precision'] > 0.9


@pytest.mark.parametrize('working_directory', ['ibis-king'], indirect=True)
def test_simulation_ibis_king(docker_client, grape_image, reference_directory, working_directory):
    volumes = {
        reference_directory.path: {'bind': CONTAINER_REFERENCE_DIRECTORY, 'mode': 'ro'},
        working_directory: {'bind': CONTAINER_WORKING_DIRECTORY, 'mode': 'rw'}
    }

    simulate_command = f'launcher.py simulate --ref-directory {CONTAINER_REFERENCE_DIRECTORY} --cores 8 ' \
                       f'--directory {CONTAINER_WORKING_DIRECTORY} --flow ibis-king --assembly hg37 ' \
                       f'--seed 42 --real-run'

    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=simulate_command, volumes=volumes)

    # Read file result with simulation metrics
    metrics = _read_metrics_file(os.path.join(working_directory, METRICS_FILEPATH))

    # Validate simultation metrics
    assert metrics['1']['Recall'] > 0.99 and metrics['1']['Precision'] > 0.99
    assert metrics['2']['Recall'] > 0.99 and metrics['2']['Precision'] > 0.99
    assert metrics['3']['Recall'] > 0.99 and metrics['3']['Precision'] > 0.97

    assert metrics['4']['Recall'] > 0.90 and metrics['4']['Precision'] > 0.95
    assert metrics['5']['Recall'] > 0.90 and metrics['5']['Precision'] > 0.94
    assert metrics['6']['Recall'] > 0.80 and metrics['6']['Precision'] > 0.90

    assert metrics['7']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['8']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['9']['Recall'] > 0 and metrics['1']['Precision'] > 0.9


@pytest.mark.parametrize('working_directory', ['germline-king'], indirect=True)
def test_simulation_germline_king(docker_client, grape_image, reference_directory, working_directory):
    volumes = {
        reference_directory.path: {'bind': CONTAINER_REFERENCE_DIRECTORY, 'mode': 'ro'},
        working_directory: {'bind': CONTAINER_WORKING_DIRECTORY, 'mode': 'rw'}
    }

    simulate_command = f'launcher.py simulate --ref-directory {CONTAINER_REFERENCE_DIRECTORY} --cores 8 ' \
                       f'--directory {CONTAINER_WORKING_DIRECTORY} --flow germline-king ' \
                       f'--phase --impute --assembly hg37 --real-run --seed 42 '

    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=simulate_command, volumes=volumes)

    # Read file result with simulation metrics
    metrics = _read_metrics_file(os.path.join(working_directory, METRICS_FILEPATH))

    # Validate simultation metrics
    assert metrics['1']['Recall'] > 0.99 and metrics['1']['Precision'] > 0.99
    assert metrics['2']['Recall'] > 0.99 and metrics['2']['Precision'] > 0.99
    assert metrics['3']['Recall'] > 0.99 and metrics['3']['Precision'] > 0.98

    assert metrics['4']['Recall'] > 0.90 and metrics['4']['Precision'] > 0.95
    assert metrics['5']['Recall'] > 0.90 and metrics['5']['Precision'] > 0.95
    assert metrics['6']['Recall'] > 0.79 and metrics['6']['Precision'] > 0.90

    assert metrics['7']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['8']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['9']['Recall'] > 0 and metrics['1']['Precision'] > 0.9


@pytest.mark.parametrize('working_directory', ['khazar'], indirect=True)
def test_khazar(docker_client, grape_image, reference_directory, working_directory, test_data_directory):
    volumes = {
        reference_directory.path: {'bind': CONTAINER_REFERENCE_DIRECTORY, 'mode': 'ro'},
        working_directory: {'bind': CONTAINER_WORKING_DIRECTORY, 'mode': 'rw'},
        test_data_directory: {'bind': CONTAINER_TESTING_REAL_DATA_DIRECTORY, 'mode': 'ro'}
    }

    preprocess_command = f'launcher.py preprocess --ref-directory {CONTAINER_REFERENCE_DIRECTORY} --cores 8 ' \
                         f'--directory {CONTAINER_WORKING_DIRECTORY} --vcf-file {KHAZAR_VCF} ' \
                         f'--assembly hg37 --real-run'

    find_command = f'launcher.py find --ref-directory {CONTAINER_REFERENCE_DIRECTORY} --cores 8 ' \
                   f'--directory {CONTAINER_WORKING_DIRECTORY} --flow ibis --real-run'

    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=preprocess_command, volumes=volumes)
    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=find_command, volumes=volumes)

    # Read file result with relatives
    relatives = _read_relatives_file(os.path.join(working_directory, RELATIVES_FILEPATH))

    # Validate relatives
    num_of_relatives = len(list(relatives))
    assert 55 <= num_of_relatives <= 57


@pytest.mark.parametrize('working_directory', ['aadr'], indirect=True)
def test_aadr(docker_client, grape_image, reference_directory, working_directory, test_data_directory):
    volumes = {
        reference_directory.path: {'bind': CONTAINER_REFERENCE_DIRECTORY, 'mode': 'ro'},
        working_directory: {'bind': CONTAINER_WORKING_DIRECTORY, 'mode': 'rw'},
        test_data_directory: {'bind': CONTAINER_TESTING_REAL_DATA_DIRECTORY, 'mode': 'ro'}
    }

    preprocess_command = f'launcher.py preprocess --ref-directory {CONTAINER_REFERENCE_DIRECTORY} --cores 8 ' \
                         f'--directory {CONTAINER_WORKING_DIRECTORY} --vcf-file {AADR_VCF} --het-samples 0.0 ' \
                         f'--assembly hg37 --real-run'

    find_command = f'launcher.py find --ref-directory {CONTAINER_REFERENCE_DIRECTORY} --cores 8 ' \
                   f'--directory {CONTAINER_WORKING_DIRECTORY} --flow ibis --real-run'

    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=preprocess_command, volumes=volumes)
    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=find_command, volumes=volumes)

    # Read file result with relatives
    relatives = _read_relatives_file(os.path.join(working_directory, RELATIVES_FILEPATH))

    # Validate relatives
    aadr_samples = _read_samples_file(AADR_SAMPLES_FILEPATH)
    failed_relations = 0
    for relationship in relatives:
        if aadr_samples[relationship[0]] != aadr_samples[relationship[1]]:
            failed_relations += 1
    assert failed_relations < 2
    
