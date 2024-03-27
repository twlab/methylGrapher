#!/usr/bin/env python

import sys
import distutils

# Get version automatically
sys.path.append('./src')
import main

ver = main.__version__

distutils.core.setup(
      name='methylGrapher',
      version=ver,
      description='WGBS data analysis for genome graph',
      author='Wenjin Zhang',
      author_email='wenjin@wustl.edu',
      url='https://github.com/twlab/methylGrapher',
      download_url=f'https://github.com/twlab/methylGrapher/releases/download/V{ver}/methylGrapher-{ver}.tar.gz',

      python_requires='>=3.7',

      scripts=[
            'bin/methylGrapher',
      ],

      install_requires=[
            'pgzip',
      ],

      classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Developers',
            'Topic :: Software Development :: Build Tools',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10'
      ],

)





