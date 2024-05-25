# Supervised Learning for Biomarker Classification in OSF
This section explores the application of supervised learning for classifying novel genes into predefined biomarker classes within the context of OSF scRNA-seq data. The goal is to develop a model that can predict whether a newly discovered gene exhibits overexpressed, underexpressed, or normal expression patterns within specific cell types.

## Data Preparation
* **Feature Selection**: The gene expression data (obtained from scRNA-seq) serves as the primary feature set for model training.
* Cell Type Information: Cell type labels, assigned through UMAP visualization and cell type annotation, are incorporated as additional features for some models.
* Data Integration: The gene expression data and cell type information are combined into a single dataset for model training.
* Preprocessing: The data may undergo normalization and scaling procedures to ensure compatibility with the chosen machine learning algorithms.

## Models Explored
This section evaluates the performance of various supervised learning algorithms for the classification task:
* **XGBoost (Extreme Gradient Boosting)**: A powerful ensemble learning method known for its speed and accuracy.
* **LightGBM (Light Gradient Boosting Machine)**: A highly efficient algorithm utilizing gradient-based learning and feature importance weighting.
* **K-Nearest Neighbors (KNN)**: A simple and interpretable algorithm that classifies data points based on their nearest neighbors in the training set.
* **Mixture of Experts (MoE)**: An advanced technique combining multiple "expert" models (e.g., XGBoost, LightGBM, KNN) into a single model, leveraging the strengths of each for improved performance.
* **Custom Biomarker Classification Model**: We propose a novel model that integrates biological knowledge gleaned from scRNA-seq data. This model employs a two-level classification approach based on observed functional relationships between cell type clusters.

## Model Training and Evaluation
* Each model is trained using a portion of the prepared data.
* Accuracy serves as the primary metric for evaluating model performance.
* Additional metrics may be employed depending on the specific model and classification task.

## Results
This section presents the classification accuracy achieved by each model, along with any relevant visualizations or comparisons. The performance of our custom model incorporating biological insights will be compared to the established machine learning algorithms. We evaluated the performance of each model using accuracy as the primary metric. Here's a summary of the results:

* *XGBoost*: Achieved an accuracy of 47.42%, but could only identify 5 out of 8 clusters, limiting its applicability.

* *LightGBM*: Offered improved cluster identification (all eight clusters) with an accuracy of 46.85%.

* *KNN*: Exhibited the lowest accuracy 43.6%, likely due to its reliance on nearest neighbors rather than gradient descent-based learning. However, it also identified all eight clusters.

* *MoE*: Achieved an accuracy of 47.83% and successfully identified all eight clusters. The gating network likely played a crucial role in leveraging the strengths of each individual expert model.

* *Custom Biomarker Classification Model*: Notably, our custom model achieved a competitive accuracy of 53.94% and successfully identified all eight clusters.  This performance suggests that incorporating biological knowledge into the model development process can be beneficial.

## Future Directions
This section discusses potential avenues for further exploration and refinement:
* Hyperparameter optimization for each model.
* Feature engineering techniques to improve classification performance.
* Integration of additional domain knowledge into the model development process.

This README provides a high-level overview of the supervised learning section in your research project. Remember to replace the placeholders with specific details from your analysis.  Feel free to expand upon each section with more specific information relevant to your implementation.