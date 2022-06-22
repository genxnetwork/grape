import os
import pytest
import docker
import csv
import shutil

from datetime import datetime
from reference_directory import ReferenceDirectory


HOME_DIRECTORY = os.path.expanduser('/media/')

GRAPE_DOCKERFILE = 'containers/snakemake/Dockerfile'
GRAPE_IMAGE_TAG = 'genx_relatives:latest'

REFERENCE_DIRECTORY = os.path.join(HOME_DIRECTORY, 'ref_test')
CONTAINER_REFERENCE_DIRECTORY = '/media/ref'
CONTAINER_WORKING_DIRECTORY = '/media/data'
TESTING_REAL_DATA_DIRECTORY = '/media/test_data'
KHAZAR_VCF = os.path.join(TESTING_REAL_DATA_DIRECTORY, 'khazar314.vcf.gz')
AADR_VCF = os.path.join(TESTING_REAL_DATA_DIRECTORY, 'aadr.reheaded.vcf.gz')
METRICS_FILEPATH = 'results/metrics.tsv'
RELATIVES_FILEPATH = 'results/relatives.tsv'
AADR_SAMPLES_FILEPATH = os.path.join(TESTING_REAL_DATA_DIRECTORY, 'aadr_samples.csv')


def _get_download_reference_command(reference_directory):
    return f'launcher.py reference --use-bundle --ref-directory {reference_directory} ' \
            '--phase --impute --real-run'


def _get_simulate_command(reference_directory, working_directory, flow):
    return f'launcher.py simulate --ref-directory {reference_directory} --cores 8 ' \
           f'--directory {working_directory} --flow {flow} --assembly hg37 --seed 42 --real-run'


def _get_preprocess_command(reference_directory, working_directory, input_file):
    return f'launcher.py preprocess --ref-directory {reference_directory} --cores 8 ' \
           f'--directory {working_directory} --vcf-file {input_file} --assembly hg37 --real-run'


def _get_find_command(reference_directory, working_directory):
    return f'launcher.py find --ref-directory {reference_directory} --cores 8 ' \
           f'--directory {working_directory} --flow ibis --real-run'


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
        reader = csv.DictReader(relatives_file, delimiter="\t")
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

    yield docker_client.images.get(GRAPE_IMAGE_TAG)

    # Fixture teardown to remove GRAPE Docker image
    docker_client.images.remove(GRAPE_IMAGE_TAG, force=True, noprune=False)


@pytest.fixture
def reference_directory(docker_client, grape_image) -> ReferenceDirectory:
    reference_directory = ReferenceDirectory(REFERENCE_DIRECTORY)

    if not reference_directory.is_valid():
        command = _get_download_reference_command(CONTAINER_REFERENCE_DIRECTORY)
        volumes = {
            reference_directory.path: {'bind': CONTAINER_REFERENCE_DIRECTORY, 'mode': 'rw'}
        }

        docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=command, volumes=volumes)

    return reference_directory


@pytest.fixture
def test_data_directory():
    return TESTING_REAL_DATA_DIRECTORY


@pytest.fixture(scope="function")
def working_directory(request):
    utc_timestamp = datetime.utcnow().strftime("%Y%m%d-%H%M%S-utc")
    working_directory_name = '-'.join([request.param, utc_timestamp])
    working_directory_path = os.path.join(HOME_DIRECTORY, working_directory_name)

    yield working_directory_path

    # Fixture teardown to remove working directory
    shutil.rmtree(working_directory_path)


@pytest.fixture(scope="function")
def simulate_command(request):
    return _get_simulate_command(CONTAINER_REFERENCE_DIRECTORY, CONTAINER_WORKING_DIRECTORY, request.param);


@pytest.fixture()
def find_command():
    return _get_find_command(CONTAINER_REFERENCE_DIRECTORY, CONTAINER_WORKING_DIRECTORY);


@pytest.fixture(scope="function")
def preprocess_command(request):
    return _get_preprocess_command(CONTAINER_REFERENCE_DIRECTORY, CONTAINER_WORKING_DIRECTORY, request.param);


simulation_list = [('ibis', 'simulation-ibis'),
                   ('ibis-king', 'simulation-ibis-king'),
                   ('germline-king --assembly hg38 --phase --impute',  # germline needs extra flags to work
                    'simulation-germline-king')]

real_data_list = [('simulation-khazar', KHAZAR_VCF),
                  ('simulation-aadr', f'{AADR_VCF} --het-samples 0.0')]  # ancient samples have zero heterozygosity


@pytest.mark.parametrize('simulate_command,working_directory', simulation_list, indirect=True)
def test_simulation(docker_client, grape_image, reference_directory, working_directory, simulate_command):
    volumes = {
        reference_directory.path: {'bind': CONTAINER_REFERENCE_DIRECTORY, 'mode': 'ro'},
        working_directory: {'bind': CONTAINER_WORKING_DIRECTORY, 'mode': 'rw'}
    }

    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=simulate_command, volumes=volumes)

    # Read file result with simulation metrics
    metrics = _read_metrics_file(os.path.join(working_directory, METRICS_FILEPATH))

    # Validate simultation metrics
    assert metrics['1']['Recall'] > 0.99 and metrics['1']['Precision'] > 0.99
    assert metrics['2']['Recall'] > 0.99 and metrics['2']['Precision'] > 0.99
    assert metrics['3']['Recall'] > 0.99 and metrics['3']['Precision'] > 0.98

    assert metrics['4']['Recall'] > 0.90 and metrics['4']['Precision'] > 0.95
    assert metrics['5']['Recall'] > 0.90 and metrics['5']['Precision'] > 0.95
    assert metrics['6']['Recall'] > 0.80 and metrics['6']['Precision'] > 0.90

    assert metrics['7']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['8']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['9']['Recall'] > 0 and metrics['1']['Precision'] > 0.9


@pytest.mark.parametrize('working_directory,preprocess_command', real_data_list, indirect=True)
def test_real_data(docker_client, grape_image, reference_directory,
                   working_directory, test_data_directory, find_command, preprocess_command):
    volumes = {
        reference_directory.path: {'bind': CONTAINER_REFERENCE_DIRECTORY, 'mode': 'ro'},
        working_directory: {'bind': CONTAINER_WORKING_DIRECTORY, 'mode': 'rw'},
        test_data_directory: {'bind': TESTING_REAL_DATA_DIRECTORY, 'mode': 'ro'}
    }

    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=preprocess_command, volumes=volumes)
    docker_client.containers.run(GRAPE_IMAGE_TAG, remove=True, command=find_command, volumes=volumes)

    # Read file result with relatives
    relatives = _read_relatives_file(os.path.join(working_directory, RELATIVES_FILEPATH))

    # Validate relatives
    current_data = working_directory.split('-')[0]
    if current_data == 'khazar':
        num_of_relatives = len(list(relatives))
        assert 55 >= num_of_relatives >= 57
    elif current_data == 'aadr':
        aadr_samples = _read_samples_file(AADR_SAMPLES_FILEPATH)
        for relative in relatives.iterrows():
            assert aadr_samples[relative[0]] == aadr_samples[relative[1]]
