"""Variant Calling Pipeline"""
import os
import subprocess
import config
import sys


def main(filenames):
    """Set dir paths to input and output files on local machine"""
    FilePath.filenames = filenames
    FilePath.set_directories()
    FilePath.set_file_paths()

    """Download reference and sequencing files from GCloud Bucket"""
    GCloud.seq_files = [filenames['seq_file_1'], filenames['seq_file_2']]
    GCloud.ref_file = filenames['ref_file']
    GCloud.ssh_login()
    GCloud.download_seq_files()
    GCloud.download_ref_genome()

    """Start NGS pipeline"""
    ngs_pipeline = NGSPipeline(filepath=FilePath)
    ngs_pipeline.bwa_index()
    ngs_pipeline.move_files()
    ngs_pipeline.bwa_mem_align()
    ngs_pipeline.samtools_sort()
    ngs_pipeline.bcftools_variant_call()
    ngs_pipeline.check_quality()


class NGSPipeline:
    """NGS pipeline written in python using following shell tools:
    bwa mem to align forward and reverse seq reads to reference genome hg19.fa
    samtools to sort aligned file .sam and convert .sam to .bam compressed binary file format
    bcftools to call variants.
    """
    def __init__(self, filepath: object):
        self.filepath = filepath
        print(self.filepath.filenames)

    def bwa_index(self):
        """Indexing the reference genome is required before aligning sequencing files"""
        try:
            subprocess.run(f"bwa index {os.path.join(config.DATA_DIR, self.filepath.filenames['ref_file'])}", stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
        except subprocess.CalledProcessError as e:
            print(e.stderr)

    def move_files(self):
        """Move index files to subfolder inside main dataset directory"""
        try:
            subprocess.run(f"mv {os.path.join(config.DATA_DIR, '*fa.*')} {os.path.join(config.DATA_DIR, 'bwa')}",
                           stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, shell=True)
        except subprocess.CalledProcessError as e:
            print(e.stderr)
        else:
            print("Files moved to subdirectory")

    def bwa_mem_align(self):
        """Align sequencing files with forward and reverse reads to reference genome.
         Outputs ASCII encoded text file of type .sam
         """
        try:
            with open(self.filepath.sam, 'w+') as output_f:
                subprocess.run(['bwa', 'mem', self.filepath.ref_file, self.filepath.seq_file_1, self.filepath.seq_file_2],
                               stdout=output_f)
        except IOError as err:
            raise err

    def samtools_sort(self):
        """Sort aligned sam file using samtools. Outputs compressed binary file."""
        try:
            with open(self.filepath.sorted_bam, 'w+') as output_f:
                subprocess.run(['samtools', 'sort', self.filepath.sam], stdout=output_f, stderr=subprocess.PIPE,
                               check=True)
        except IOError as e:
            raise e

    def bcftools_variant_call(self):
        """Call variants with bcftools"""
        cmd = f"bcftools mpileup -Ou -f {os.path.join(config.DATA_DIR, self.filepath.filenames['ref_file'])} " \
              f"{self.filepath.sorted_bam} | bcftools call -vmO v -o dataset/result.vcf"
        subprocess.run(cmd, shell=True)

    def check_quality(self):
        """Check number of variants found as a quality matrix"""
        cmd = "grep -v '##' dataset/result.vcf | wc -l"
        p = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True, text=True)
        print('Number of variants found: ', p.stdout)


class GCloud:
    """Connect to Google Cloud and download reference genome and seq files from bucket"""
    seq_files: list = None
    ref_file: str = None

    @classmethod
    def auth(cls):
        """Login to google cloud. This step is optional if --authorize-session in ss_login is set"""
        try:
            subprocess.run(['gcloud', 'auth', 'login'], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            print(e.stderr)
            sys.exit()

    @classmethod
    def ssh_login(cls):
        """Login to google cloud shell"""
        try:
            subprocess.run('gcloud cloud-shell ssh --authorize-session',
                           stdin=subprocess.PIPE,
                           check=True,
                           shell=True
                           )
        except subprocess.CalledProcessError as e:
            print(e.stderr)
            sys.exit()
        else:
            print('Gcloud ssh login success')

    @classmethod
    def download_seq_files(cls):
        for seq_file in cls.seq_files:
            print(seq_file)
            try:
                subprocess.run(
                    f"gcloud cloud-shell scp cloudshell:{os.path.join(config.GCLOUD_BUCKET, seq_file)} "
                    f"localhost:{config.DATA_DIR}", stderr=subprocess.PIPE, check=True, shell=True)
            except subprocess.CalledProcessError as e:
                print(e.stderr)

    @classmethod
    def download_ref_genome(cls):
        try:
            subprocess.run(
                f"gcloud cloud-shell scp cloudshell:~/tiny-test-data/genomes/Hsapiens/hg19/seq/{cls.ref_file} "
                f"localhost:{config.DATA_DIR}", stderr=subprocess.PIPE, check=True, shell=True)
        except subprocess.CalledProcessError as e:
            print(e.stderr)


class FilePath:
    filenames = None
    ref_file = None
    seq_file_1 = None
    seq_file_2 = None
    sam = None
    sorted_bam = None
    vcf = None

    @classmethod
    def set_directories(cls):
        """Setup folders to store downloaded sequencing/reference and pipeline output files"""
        try:
            print("Set directories")
            os.makedirs(config.BWA_DIR)
        except IOError as e:
            print("Directories already exist")
            pass
        else:
            print('Directories set')

    @classmethod
    def set_file_paths(cls):
        try:
            cls.ref_file = os.path.join(config.BWA_DIR, cls.filenames['ref_file'])
            cls.seq_file_1 = os.path.join(config.DATA_DIR, cls.filenames['seq_file_1'])
            cls.seq_file_2 = os.path.join(config.DATA_DIR, cls.filenames['seq_file_2'])
            cls.sam = os.path.join(config.DATA_DIR, cls.filenames['alignment_file'])
            cls.sorted_bam = os.path.join(config.DATA_DIR, "sorted-" + 'alignment.bam')
            cls.vcf = os.path.join(config.DATA_DIR, 'result.vcf')
        except IOError as e:
            raise e
        else:
            print("All paths set")


def set_user_input():
    ref_file = str(input('Enter name of ref_file in bucket: ') or 'hg19.fa')
    seq_file_1 = str(input('Enter name of forward sequencing file in bucket: ') or 'mt_1.fq')
    seq_file_2 = str(input('Enter name of reverse sequencing file in bucket: ') or 'mt_2.fq.gz')
    alignment_file = str(input('Enter name of sam file, output.sam: ') or 'output.sam')
    return dict(ref_file=ref_file, seq_file_1=seq_file_1, seq_file_2=seq_file_2, alignment_file=alignment_file)


if __name__ == "__main__":
    main(set_user_input())
    sys.exit()