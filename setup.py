import setuptools

setuptools.setup(
    name="MagmaPandas",
    version="1.0",
    description="...",
    author="Thomas van Gerve",
    packages=setuptools.find_packages(where="src", exclude=[]),
    package_dir={"": "src"},
    package_data={"elements": ["data/*"], "geoplot": ["data/*"]},
    install_requires=[
        "pandas",
        "scipy",
        "numpy",
        "alive_progress",
        "matplotlib-label-lines",
        "statsmodels",
    ],
)
