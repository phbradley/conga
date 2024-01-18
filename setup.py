from setuptools import setup #, find_packages
#from codecs import open
#from os import path

#here = path.abspath(path.dirname(__file__))

# Thanks to Neal Smith for getting me (Phil) started with this!

setup(
    name='conga',

    version='0.1.2',

    author='Phil Bradley',
    author_email= 'pbradley@fredhutch.org',

    license='MIT',

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English',
    ],

    install_requires=["scanpy", "leidenalg", "natsort"],

    packages=['conga','conga.tcrdist'],
    #packages=find_packages(),

    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #},
)
