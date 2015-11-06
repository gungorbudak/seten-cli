from setuptools import setup, find_packages


EXCLUDE_FROM_PACKAGES = []

setup(
    name='gsea_rbp',
    version='0.0.1',
    description='Inferring associations between RNA-binding proteins (RBPs) and gene sets based on binding signals.',
    url='https://bitbucket.org/gungorbudak/gsea_rbp/',
    author='Gungor Budak',
    author_email='gbudak@iupui.edu',
    license='MIT',
    packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'gsea_rbp = gsea_rbp.scripts.gsea_rbp:main',
            ]
        }
    )
