import setuptools

setuptools.setup(
    name= 'MagmaPandas',
    version= '0.1',
    description= '...',
    author= 'Thomas van Gerve',
    
    packages= setuptools.find_packages(
        where='src',
        exclude=[]
        ),
    package_dir={"": "src"},

    package_data= {"elements": ['data/*']}, 

    install_requires= [
    'pandas'
    ]
)