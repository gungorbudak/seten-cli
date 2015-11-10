from setuptools import setup, find_packages


REQUIRES = [
    'argparse>=1.2.1',
    'requests>=2.8.1',
    'intervaltree>=2.1.0',
    'numpy>=1.9.2',
    'scipy>=0.7.2'
    ]

EXCLUDE_FROM_PACKAGES = []

setup(
    name='seten',
    version='0.0.1',
    description='Gene set enrichment on CHIP-seq RBA-binding protein binding signals datasets',
    url='https://bitbucket.org/gungorbudak/seten',
    author='Gungor Budak',
    author_email='gbudak@iupui.edu',
    license='MIT',
    install_requires=REQUIRES,
    packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'seten = seten.cli:main',
            ]
        }
    )
