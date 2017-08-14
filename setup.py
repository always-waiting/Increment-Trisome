from setuptools import setup

setup(
    name="berry_increment_trisome",
    version="0.0.1",
    description="The Method for xromate trisome detection",
    author="bixichao",
    author_email="bixichao@berrygenomics.com",
    license='MIT',
    packages=['berry_increment_trisome'],
    install_requires=[
        "pandas",
        "numpy",
        "pymongo"
    ],
    include_package_data=True,
    scripts=[
        "bin/berry-increment-trisome-auto",
        "bin/berry-increment-trisome-import"
    ],
    zip_safe=False
)
