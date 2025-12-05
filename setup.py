from setuptools import setup, find_packages

setup(
    name="sapling-phylogeny",  # PyPI name (sapling is taken)
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "tqdm",
    ],
    entry_points={
        "console_scripts": [
            # This creates the 'sapling' command in the terminal
            # syntax: command_name = package.file:function
            "sapling=sapling.sapling:main",
        ],
    },
)
