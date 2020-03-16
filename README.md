
# MAP fingerprint - Design and Documentation  

The canonical, not isomeric, and rooted SMILES of the circular substructures `CS` from radius one up to a user-given radius `n` (default `n=2`, `MAP4`) are generated for each atom. All atom pairs are extracted, and their minimum topological distance `TP` is calculated. For each atom pair `jk`, for each considered radius `r`, a `Shingle` is encoded as: `CS`<sub>`rj`</sub>`|TP`<sub>`jk`</sub>`|CS`<sub>`rk`</sub> , where the two `CS` are annotated in alphabetical order, resulting in n Shingles for each atom pairs. 

![MAP4 atom pair encoding scheme](https://cloud.gdb.tools/s/oANAxRazApL5EDw/preview)

The resulting list of Shingles is hashed using the unique mapping `SHA-1` to a set of integers `S`<sub>`i`</sub>, and its correspondent transposed vector `s`<sup>`T`</sup><sub>`i`</sub> is MinHashed.

![MihHash](https://cloud.gdb.tools/s/nLjQKTcHPLdpnxJ/preview)

To use the MAP4 fingerprint:
- `git clone <repo url>`
- `cd <repo folder>`

A Conda environment.yml is supplied with all the required libraries:
- `conda env create -f environment.yml`
- `conda activate map4`

Run the fingerprint from terminal
- `python map4.py -i smilesfile.smi -o outputfile`

Or import the MAP4Calculator class in your python file (see `test.py`)

<pre>


</pre>
# MAP4 - Similarity Search of ChEMBL, Human Metabolome, and SwissProt

Draw a structure or paste its SMILES, or write a natural peptides linear sequence.
Search for its analogs in the MAP4 or MHFP6 space of ChEMBL, of the Human Metabolome Database (HMDB), or of the 'below 50 residues subset' of SwissProt.

The MAP4 search can be found at: http://map-search.gdb.tools/.

The code of the MAP4 similarity search can be found in this repository folder `MAP4-Similarity-Search`

To run the app locally:
- Download the MAP4SearchData [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3671214.svg)](https://doi.org/10.5281/zenodo.3671214)
- Run `docker run -p 8080:5000 --mount type=bind,target=/MAP4SearchData,source=/your/absolut/path/MAP4SearchData  --restart always --name mapsearch alicecapecchi/map-search:latest`
- The MAP4 similarity search will be running at http://0.0.0.0:8080/

<pre>


</pre>

# Extended Benchmark

Compounds and training list used to extend the Riniker et. al. fingerprint benchmark (Riniker, G. Landrum, J. Cheminf., 5, 26 (2013), DOI: 10.1186/1758-2946-5-26, URL: http://www.jcheminf.com/content/5/1/26, GitHub page: https://github.com/rdkit/benchmarking_platform) to peptides.
