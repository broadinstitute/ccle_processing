from depmapomics.qc.config import VIRTUAL_RELEASE, REFERENCE_RELEASE

def pytest_sessionstart(session):
    print('REFERENCE_RELEASE: {}.{}'.format(REFERENCE_RELEASE['name'], REFERENCE_RELEASE['version']))
    print('VIRTUAL_RELEASE: {}.{}'.format(VIRTUAL_RELEASE['name'], VIRTUAL_RELEASE['version']))
