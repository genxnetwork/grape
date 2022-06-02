import os
import pytest
import docker
import csv
import shutil

from datetime import datetime
from reference_directory import ReferenceDirectory


HOME_DIRECTORY = os.path.expanduser('~')

GRAPE_DOCKERFILE = 'containers/snakemake/Dockerfile'
GRAPE_IMAGE_TAG = 'genx_relatives:latest'

REFERENCE_DIRECTORY = os.path.join(HOME_DIRECTORY, 'ref')
CONTAINER_REFERENCE_DIRECTORY = '/media/ref'
CONTAINER_WORKING_DIRECTORY = '/media/data'
METRICS_FILEPATH = 'results/metrics.tsv'


def _get_download_reference_command(reference_directory):
    return f'launcher.py reference --use-bundle --ref-directory {reference_directory} ' \
            '--phase --impute --real-run'


def _get_simulate_command(reference_directory, working_directory):
    return f'launcher.py simulate --ref-directory {reference_directory} --cores 8 ' \
           f'--directory {working_directory} --flow ibis --assembly hg37 --seed 42 --real-run'


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
def working_directory():
    utc_timestamp = datetime.utcnow().strftime("%Y%m%d-%H%M%S-utc")
    working_directory_name = '-'.join('simultation-ibis', utc_timestamp)
    working_directory_path = os.path.join(HOME_DIRECTORY, working_directory_name)

    yield working_directory_path

    # Fixture teardown to remove working directory
    shutil.rmtree(working_directory_path)


@pytest.fixture
def simulate_command():
    return _get_simulate_command(CONTAINER_REFERENCE_DIRECTORY, CONTAINER_WORKING_DIRECTORY);


def test_simulation(docker_client, grape_image, reference_directory, working_directory, simulate_command):
    working_directory = os.path.join(HOME_DIRECTORY, 'simulation-ibis')
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
    assert metrics['3']['Recall'] > 0.99 and metrics['3']['Precision'] > 0.99

    assert metrics['4']['Recall'] > 0.90 and metrics['4']['Precision'] > 0.95
    assert metrics['5']['Recall'] > 0.90 and metrics['5']['Precision'] > 0.95
    assert metrics['6']['Recall'] > 0.80 and metrics['6']['Precision'] > 0.90

    assert metrics['7']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['8']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
    assert metrics['9']['Recall'] > 0 and metrics['1']['Precision'] > 0.9
