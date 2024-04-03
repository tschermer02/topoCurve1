from setuptools import setup, find_packages

setup(
    name='topo_curve',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'Pillow',
        'photutils',
        'tifffile',
        'geotiff', 
        'scikit-learn'
    ],
    python_requires='>=3.6',
    author='Taylor Schermer, Joel Nash, Nate Klema',
    author_email='tschermer@fortlewis.edu, jxnash@fortlewis.edu, ntklema@fortlewis.edu',
    description='A package for processing digital elevation models.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/tschermer02/topoCurve1',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
