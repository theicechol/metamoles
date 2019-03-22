from setuptools import setup, find_packages

setup(
    name = 'metamoles',
    version = '1.0.0',
    url = 'https://github.com/theicechol/metamoles.git',
    author = 'theicechol',
    author_email = 'cholpisit@gmail.com’,
    description = 'UWDIRECT-MetaMolES',
    packages = find_packages(),
    install_requires = ['numpy >= 1.11.1', 'matplotlib >= 1.5.1'],
)
