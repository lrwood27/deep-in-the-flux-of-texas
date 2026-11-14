# Input and processed data

- `raw/boundaries/`: original Census boundary ZIP, tracked with its embedded metadata.
- `raw/population/`: supplied population-density GeoTIFF, tracked; original attribution remains to be recorded.
- `raw/elevation/urls.txt`: tracked list of USGS tile URLs.
- `raw/elevation/tiles/`: local downloads and partial files, ignored by Git.
- `processed/`: generated mesh/region arrays for solver input.

Keep original inputs intact. Scripts resolve these locations from their own file paths.
Future cross-section datasets should have their energy, units, provenance, density,
and composition documented before being used; their format/source is undecided.
