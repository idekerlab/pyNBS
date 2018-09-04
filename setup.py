"""
Setup module adapted from setuptools code. See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

setup(
	name='pyNBS',
	version='0.2.0',
	description='Python package to perform network based stratification of binary somatic mutations as described in Hofree et al 2013.',
	url='https://github.com/huangger/pyNBS',
	author='Justin Huang',
	author_email='jkh013@ucsd.edu',
	license='MIT',
	classifiers=[
		'Development Status :: 4 - Beta',
		'Intended Audience :: Science/Research',
		'Topic :: Software Development :: Build Tools',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2.7'
	],
	packages=find_packages(exclude=['os', 'random', 'time']),
	install_requires=[
        'lifelines>=0.9.1',
        'networkx>=2.0',
        'numpy>=1.11.0',
        'matplotlib>=1.5.1',
        'pandas>=0.19.0',
        'scipy>=0.17.0',
        'scikit-learn>=0.17.1',
        'seaborn>=0.7.1']
)
