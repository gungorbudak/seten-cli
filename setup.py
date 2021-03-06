from setuptools import setup, find_packages


NUMPY = 'numpy>=1.11.0'

REQUIRES = [
    'argparse>=1.2.1',
    'requests>=2.8.1',
    'intervaltree>=2.1.0',
    'scipy>=0.17.0',
    NUMPY
]

EXCLUDE_FROM_PACKAGES = []

setup(
    name='seten',
    version='0.0.3',
    description='Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets',
    url='https://github.com/gungorbudak/seten-cli',
    author='Gungor Budak',
    author_email='gbudak@iupui.edu',
    license='MIT',
    setup_requires=[NUMPY],
    install_requires=REQUIRES,
    packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'seten = seten.cli:main',
        ]
    }
)
