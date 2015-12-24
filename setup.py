from setuptools import setup, find_packages


REQUIRES = [
    'argparse>=1.2.1',
    'requests>=2.8.1',
    'intervaltree>=2.1.0',
    'numpy>=1.10.1',
    'scipy>=0.16.1'
]

EXCLUDE_FROM_PACKAGES = []

setup(
    name='seten',
    version='0.0.1',
    description='Gene set enrichment on CLIP-seq RNA-binding protein binding signals datasets',
    url='https://github.com/gungorbudak/seten-cli',
    author='Gungor Budak',
    author_email='gbudak@iupui.edu',
    license='MIT',
    install_requires=REQUIRES,
    packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
    package_data={
        'seten': [
            'resources/*'
            'resources/collections/*'
            ]
        },
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'seten = seten.cli:main',
        ]
    }
)
