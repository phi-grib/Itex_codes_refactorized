from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
   name='itex_codes_refactorized',
   version='1.0',
   description='Codes to mantain CII refactorized version',
   license='GNU',
   long_description=long_description,
   author='Eric March Vila, Giacomo Ferretti',
   author_email='eric.march@upf.edu','gferretti85@gmail.com',
   url='https://github.com/phi-grib/Itex_codes_refactorized',
   packages=['UpdateDB']
)
