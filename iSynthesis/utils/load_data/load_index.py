import faiss
from pickle import load, dump

index = faiss.read_index('iSynthesis/utils/load_data/trained_cosine_from_python3.8_may.index', faiss.IO_FLAG_ONDISK_SAME_DIR)
ids = load(open('iSynthesis/utils/load_data/db_id.pickle', 'rb'))