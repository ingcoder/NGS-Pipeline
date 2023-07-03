import os

HOME_DIR = os.environ["HOME"]
PROJECT_FOLDER = "path/to/project/folder" 
PROJECT_DIR = os.path.join(HOME_DIR, PROJECT_FOLDER)
DATA_FOLDER = "dataset"
BWA_FOLDER = 'bwa'
DATA_DIR = os.path.join(PROJECT_DIR, DATA_FOLDER)
BWA_DIR = os.path.join(PROJECT_DIR, DATA_FOLDER, BWA_FOLDER)
GCLOUD_BUCKET = os.path.join('~', 'you_dataset_folder', 'subdir')


