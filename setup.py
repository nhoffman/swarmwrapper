import os
import subprocess
from setuptools import setup, find_packages

subprocess.call(
    ('mkdir -p swarmwrapper/data && '
     'git describe --tags --dirty > swarmwrapper/data/ver.tmp'
     '&& mv swarmwrapper/data/ver.tmp swarmwrapper/data/ver '
     '|| rm -f swarmwrapper/data/ver.tmp'),
    shell=True, stderr=open(os.devnull, "w"))

from swarmwrapper.swarmwrapper import __version__

setup(
    author='Noah Hoffman',
    author_email='noah.hoffman@gmail.com',
    description='Wrapper for virtualenv, pip, and wheel',
    url='https://github.com/nhoffman/swarmwrapper',
    name='swarmwrapper',
    packages=find_packages(),
    package_dir={'swarmwrapper': '.'},
    package_data={'swarmwrapper': ['data/ver']},
    entry_points={'console_scripts': ['swarmwrapper = swarmwrapper:main']},
    version=__version__,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2.7',
    ],
)
