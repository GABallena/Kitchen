import tensorflow as tf
import numpy as np

def extract_features(alignment):
    # Example: Encode each column as a one-hot matrix
    residues = "ACDEFGHIKLMNPQRSTVWY-"
    residue_to_index = {res: i for i, res in enumerate(residues)}
    columns = list(zip(*alignment))
    features = [
        [residue_to_index.get(res, len(residues)) for res in col] for col in columns
    ]
    return np.array(features)

def predict_alignment_quality(alignment, model_path="model.h5"):
    model = tf.keras.models.load_model(model_path)
    features = extract_features(alignment)
    features = np.expand_dims(features, axis=0)  # Add batch dimension
    prediction = model.predict(features)
    return prediction[0, 0]  # Assuming scalar output

# Load alignment file
def load_alignment(file_path):
    with open(file_path) as f:
        return [line.strip() for line in f if not line.startswith(">")]

alignment = load_alignment("input.fasta")
quality_score = predict_alignment_quality(alignment, "pretrained_model.h5")

with open("output_deep_learning_score.txt", "w") as f:
    f.write(f"Deep Learning Quality Score: {quality_score}\n")
