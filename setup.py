from setuptools import setup

setup(
    name='CaSim',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    packages=['casim',
             ],
    license='MIT',
    long_description=open('README.md').read(),
)
