'''
@author: Marta Luksza
'''

from setuptools import setup


setup(
        name="cfit",
        version="1.0",
        author="Marta Luksza",
        author_email="lukszam@mskcc.org",
        description=("Cancer fitness modeling"),
        license="MIT",
        keywords="Fitness model, neontigens",
        url="",

        package_data={'cfit.data.matrices': ['*.json', '*.txt']},
        include_package_data=True,
        packages=['cfit', 'cfit.data', 'cfit.data.matrices',
                  'cfit.fitness', 'cfit.fitness.neo_quality',
                  'cfit.tree',
                  'cfit.tree.mutation', 'cfit.tree.node',
                  'cfit.patient', 'cfit.util', 'cfit.plot'],
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Science"
            ],
        scripts=[], install_requires=['pandas', 'numpy', 'matplotlib', 'lifelines', 'Bio', 'scikit-learn', 'plotly==5.14.1']
        )
