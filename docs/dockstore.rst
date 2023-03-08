Dockstore
================================

GRAPE supports execution from the Dockstore, which is a free and open source platform for sharing reusable and scalable
analytical tools and workflows. It’s developed by the Cancer Genome Collaboratory
and implements `GA4GH Tool Registry Service API (TRS) <https://www.ga4gh.org/news/tool-registry-service-api-enabling-an-interoperable-library-of-genomics-analysis-tools/>`_.

To deploy GRAPE from the Dockstore:

1. `Install Dockstore <https://dockstore.org/quick-start>`_.
2. Clone the `GRAPE repository <https://dockstore.org/organizations/GenX/collections/GRAPE>`_.
3. Download reference datasets.

There are four available options:
::

    # Download minimal reference datasets from public available sources (without phasing and imputation)
    dockstore tool launch --local-entry workflows/reference/cwl/ref_min.cwl --json workflows/reference/cwl/config.json --script

    # Download complete reference datasets from public available sources (for phasing and imputation)
    dockstore tool launch --local-entry workflows/reference/cwl/ref.cwl --json workflows/reference/cwl/config.json --script

    # Download minimal reference datasets as bundle from (without phasing and imputation)
    dockstore tool launch --local-entry workflows/bundle/cwl/bundle_min.cwl --json workflows/reference/bundle/config.json --script

    # Download complete reference datasets as bundle from (for phasing and imputation)
    dockstore tool launch --local-entry workflows/bundle/cwl/bundle.cwl --json workflows/reference/bundle/config.json --script


4. Run the preprocessing.
::

    dockstore tool launch --local-entry workflows/preprocess2/cwl/preprocess.cwl --json workflows/preprocess2/cwl/config.json --script


5. Run the relatedness inference workflow.
::

    dockstore tool launch --local-entry workflows/find/cwl/find.cwl --json workflows/find/cwl/config.json --script


Each pipeline step requires `config.json`.

Each config has predefined directories and paths, but you can override these paths by changing them in these files.

Also notice that Dockstore saves its runtime in `datastore` directory.
This directory will grow up with each run.
We recommend to clean it up after each step, especially after reference downloading.