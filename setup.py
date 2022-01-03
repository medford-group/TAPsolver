from setuptools import setup

setup(name='tapsolver',
      version='0.1',
      description='A simulation and analysis tool for TAP reactor systems',
      url='http://github.com/medford-group/TAPsolver',
      author='Adam Yonge',
      author_email='ayonge3@gatech.edu',
      license='MIT',
      packages=find_packages(),#['tapsolver'],
      install_requires=['fenics','dolfin-adjoint','imageio','pandas','matplotlib'],
      zip_safe=False)
