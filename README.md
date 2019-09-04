VISplicing
=========

Visualization toolkit for splicing QTL and circularization QTL

## Requirements:

- Python 2.7
- pandas package
- pysam package
- matplotlib, scipy, and numpy
- samtools must be installed and in your $PATH
- tabix must be installed and in your $PATH

## Examples

### step 1: generate pickle file
```bash
python initialize_data.py chr18:9826168 chr18:9930662-9951252 --vcf test_files/SNP.chr18_9826168.vcf.gz --gtf test_files/HG19.gtf --mf test_files/map_file.txt --circ test_files/chr18_9931806_9950565.circRNA.rawcount.txt --output test_files
```
### step 2: make Sashimi plot
```bash
python plot.py test_files/chr18:9826168@chr18:9930662-9951252.p pickle settings_file
```
## Sashimi plot

<p align="center">
<img src="https://github.com/jiajiepeng/VISplicing/blob/master/plots/chr18_sashimi.jpg" width="600"/>
</p>


## Credits
We would like to thank for the open-sourcing [SplicePlot](https://github.com/wueric/SplicePlot) code.
