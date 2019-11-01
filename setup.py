from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

with open('requirements.txt') as f:
    lines = f.readlines()
install_requirements = [i.strip() for i in lines]


setup(name='concavity',
      version='0.0.1',
      description='concave hull delineation and geometry concavity and convexity detection',
      long_description=readme(),
      classifiers=[
          'Development Status :: ',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.6',
          'Topic :: GIS :: Geometry :: Concave :: Convex',
      ],
      keywords='concave hull convex geometry gis',
      url='http://github.com/mlichter2/concavity',
      author='mlichter2',
      author_email='mlichter@gmail.com',
      license='MIT',
      packages=['concavity'],
      install_requires=install_requirements,
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      include_package_data=True,
      zip_safe=False)
