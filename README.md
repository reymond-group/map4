Folder description:
- `Extended-Benchmark`: compounds and query lists used for the peptide benchmark
- `MAP4-Similarity-Search`: source code for the similarity search app
- `map4`: MAP4 fingerprint source code
 

# MAP fingerprint - Design and Documentation  

The canonical, not isomeric, and rooted SMILES of the circular substructures `CS` from radius one up to a user-given radius `n` (default `n=2`, `MAP4`) are generated for each atom. All atom pairs are extracted, and their minimum topological distance `TP` is calculated. For each atom pair `jk`, for each considered radius `r`, a `Shingle` is encoded as: `CS`<sub>`rj`</sub>`|TP`<sub>`jk`</sub>`|CS`<sub>`rk`</sub> , where the two `CS` are annotated in alphabetical order, resulting in n Shingles for each atom pairs. 

![MAP4 atom pair encoding scheme](https://cloud.gdb.tools/s/oANAxRazApL5EDw/preview)

The resulting list of Shingles is hashed using the unique mapping `SHA-1` to a set of integers `S`<sub>`i`</sub>, and its correspondent transposed vector `s`<sup>`T`</sup><sub>`i`</sub> is MinHashed.

![MihHash](https://cloud.gdb.tools/s/nLjQKTcHPLdpnxJ/preview)

To use the MAP4 fingerprint:
- `git clone https://github.com/reymond-group/map4.git@v1.0`
- `cd map4`

To install map4 trough Conda:
- `conda env create -f environment.yml`
- `conda activate map4`

To install map4 trough pip:
- install RDKit and tmap:
   - https://github.com/reymond-group/tmap
   - https://github.com/rdkit/rdkit
- `pip install git+https://github.com/reymond-group/map4@v1.0`

Run the fingerprint from terminal
- `cd map4`
- `python map4.py -i smilesfile.smi -o outputfile`

Or import the MAP4Calculator class in your python file (see `test.py`)

### Please note that the similarity/dissimilarity between two MinHashed fingerprints cannot be assessed with "standard" Jaccard, Manhattan, or Cosine functions. Due to MinHashing, the order of the features matters and the distance cannot be calculated "feature-wise". There is a well written blog post that explains it: https://aksakalli.github.io/2016/03/01/jaccard-similarity-with-minhash.html. Therefore, a custom kernel/loss function needs to be implemented for machine learning applications of MAP4 (e.g. using the distance function found in the test.py script).

# MAP4 - Similarity Search of ChEMBL, Human Metabolome, and SwissProt

Draw a structure or paste its SMILES, or write a natural peptides linear sequence.
Search for its analogs in the MAP4 or MHFP6 space of ChEMBL, of the Human Metabolome Database (HMDB), or of the 'below 50 residues subset' of SwissProt.

The MAP4 search can be found at: http://map-search.gdb.tools/.

The code of the MAP4 similarity search can be found in this repository folder `MAP4-Similarity-Search`

To run the app locally:
- Download the MAP4SearchData [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3671214.svg)](https://doi.org/10.5281/zenodo.3671214)
- Run `docker run -p 8080:5000 --mount type=bind,target=/MAP4SearchData,source=/your/absolut/path/MAP4SearchData  --restart always --name mapsearch alicecapecchi/map-search:latest`
- The MAP4 similarity search will be running at http://0.0.0.0:8080/


# Extended Benchmark

Compounds and training list used to extend the Riniker et. al. fingerprint benchmark (Riniker, G. Landrum, J. Cheminf., 5, 26 (2013), DOI: 10.1186/1758-2946-5-26, URL: http://www.jcheminf.com/content/5/1/26, GitHub page: https://github.com/rdkit/benchmarking_platform) to peptides.
