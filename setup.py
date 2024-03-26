#!/usr/bin/env python


from distutils.core import setup

setup(
      name='methylGrapher',
      version='0.1.0',
      description='WGBS data analysis for genome graph',
      author='Wenjin Zhang',
      author_email='wenjin@wustl.edu',
      url='https://github.com/twlab/methylGrapher',
      download_url='https://github.com/twlab/methylGrapher/releases/download/V0.1.0/methylGrapher-0.1.0.tar.gz',

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





