from setuptools import setup, find_packages

setup(
    name="chemem",
    version="0.0.4",
    packages=find_packages(),
    install_requires=[],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'chemem = ChemEM.chemem:main',
            'chemem.test = ChemEM.chemem:test',
            'chemem.protonate = ChemEM.chemem:protonate',
            'chemem.chemem_path = ChemEM.chemem:chimeraX_path'
        ],
    },
)
