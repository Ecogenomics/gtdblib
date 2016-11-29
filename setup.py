import os


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'gtdblib', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='gtdblib',
    version=version(),
    author='Pierre Chaumeil',
    author_email='uqpchaum@uq.edu.au',
    packages=['gtdblib'],
    package_data={'gtdblib': ['VERSION']},
    # url='http://pypi.python.org/pypi/biolib/',
    license='GPL3',
    description='Package for common tasks in GTDB and GTDBtk.',
    long_description=open('README.md').read(),
    install_requires=[],
)
