# Morse–Smale Complex for Unsupervised U-Net & Graph Learning

Topologically-aided segmentation with neural networks. This project uses the
**Morse–Smale Complex (MSC)** — a topological decomposition of an image treated
as a scalar field — as the backbone for two segmentation approaches: an
**inductive graph neural network** that classifies MSC arcs, and an
**unsupervised U-Net** trained on MSC segmentations in place of hand-drawn
ground truth.

The application throughout is blood-vessel segmentation in retinal images
(DRIVE and STARE), but the target is incidental — the study is really about
whether a topological summary can do the labeling work a human annotator
usually does.

📄 **Paper:** [`UNsup_MSC_UNET_MSC_GRAPH_LEARNING.png`](UNsup_MSC_UNET_MSC_GRAPH_LEARNING.png)
(Samuel Leventhal, 2019)

![overview](UNsup_MSC_UNET_MSC_GRAPH_LEARNING.png)

---

## The two ideas

### 1. Graph learning on the MSC

Compute the MSC over an image and you get a network of **arcs** — 1-cells
tracing ridge lines between critical points, which in a retinal scan follow the
vessels. Hand-label a small subset of arcs, then learn to classify the rest.

- Each arc carries a **feature vector** (min, max, mean, Gaussian fit, standard
  deviation) drawn from the pixels beneath it.
- Features are made more robust by projecting the arcs onto a stack of filtered
  images (difference-of-Gaussians, Sobel, Laplace, structure, blur) and
  re-sampling the same statistics under each — topology-guided augmentation.
- Learning is **inductive** (train on a region, generalize to unseen arcs),
  built on **GraphSAGE**, run on the **dual** of the MSC (arcs → nodes), with an
  unsupervised loss over random-walk node pairs and CVT-sampled seeds.
- A two-layer mean-pool aggregator (16-sample neighborhoods, 32-dim embeddings)
  plus linear regression classifies the remaining arcs.

**Result:** F1 ≈ **0.79** classifying the full held-out MSC. An ablation
(`graph_learning/`) shows a plain random forest on the raw arc features reaches
**0.909**, and DeepWalk embeddings add only 0.004 on top — evidence that the
task is geometric rather than relational.

### 2. Unsupervised U-Net

Train a **U-Net** on MSC segmentations as pseudo-ground-truth, removing the need
for hand-labels. Because "the MSC" depends on **persistence** (the prominence
threshold for surviving features), the network is given *several* segmentations
per image — persistence 7, 10, 12, 15, 20, 23 — stacked as multiple targets.

**Result:** F1 ≈ **0.46** on test (**0.28** if the low-persistence
segmentations carrying the thin vessels are dropped). Pixels that flip label
across persistence values produce contradictory targets, so the network hedges
and blurs. The proposed fix — a **spatially weighted loss** that down-weights
pixels by their instability along the MSC — is described in the paper's
conclusion.

---

## Repository layout

```
.
├── UNsup_MSC_UNET_MSC_GRAPH_LEARNING.png   the paper
├── USEME.txt                               module-level notes
│
├── MSCvisus/                               C++ Morse–Smale / discrete-gradient library
│   ├── include/                            msc* headers (mesh handlers, gradient builders, MSC)
│   └── MSCTest*/                           2D/3D viewers, basins, streaming, parallel expansion
│
├── topoml/                                 Python: MSC → graph → learning
│   ├── mscnn_seg.py                        driver: compute (Geo)MSC segmentations over images
│   ├── test_gnn.py                         driver: build MSC graph, train + classify with GraphSAGE
│   ├── graph_learning_improvement_experiment.py   driver: the CAMLP/DeepWalk/RF ablation
│   ├── graph_learning/
│   │   ├── classify_feats_graphs_v2.py     the ablation (CAMLP, DeepWalk+RF, RF variants)
│   │   ├── gbssl.py                        CAMLP label propagation
│   │   └── deepwalk/                        vendored DeepWalk
│   └── topoml/
│       ├── ArcNeuronTracer.py              build MSC/GeoMSC, encode arcs in a graph structure
│       ├── MSCGNN.py                       GraphSAGE handler: select → train → classify → view
│       ├── msc_gnn.py
│       ├── graphsage/                       GraphSAGE (aggregators, minibatch, unsupervised train)
│       ├── topology/                        MSC construction + I/O (geometric_msc, read_msc, utils)
│       ├── ui/ArcSelector.py               select arcs; draw_binary_segmentation()
│       └── eval_scripts/LinearRegression.py  classify from learned embeddings
│
└── Unsup_MSC_UNet/
    └── Unsupervised_UNet.ipynb             U-Net trained on multi-persistence MSC ground truth
```

---

## Getting started

### Dependencies

The graph-learning path was developed against **TensorFlow 1.11** (GraphSAGE),
with `scipy`, `numpy`, `matplotlib`, `scikit-learn`, `scikit-image`, and
`networkx`. The U-Net notebook uses **PyTorch**. Because the pinned versions are
old, a clean virtual environment (Python 3.7-era) is recommended.

```bash
cd topoml
pip install -r requirements.txt   # or: python setup.py develop
```

The C++ MSC library under `MSCvisus/` builds with CMake and is only needed if
you want to recompute complexes from scratch; the Python drivers can read
pre-computed MSC files.

### Typical workflow

Run the following from within `topoml/`.

**1 — Compute MSC segmentations over your images.**
Edit the data paths and `persistence_values` near the bottom of `mscnn_seg.py`,
then:

```bash
python3 mscnn_seg.py
```

This writes raw images and MSC (or Geometric MSC) segmentations per image, per
persistence.

**2 — Inductive graph learning (train + classify the MSC).**
`test_gnn.py` exposes the full pipeline — arc selection / CVT sampling,
unsupervised GraphSAGE training, classification of the remaining arcs, and a
view of the labeling:

```bash
python3 test_gnn.py
```

Key knobs live at the top of the script: `aggregator` (default
`graphsage_meanpool`), `learning_rate`, `epochs`, `depth`, the CVT sampling
tuple `(n_samples, hops)`, and the random-walk file used for the unsupervised
loss.

**3 — The graph-learning ablation** (reproduces the F1 table):

```bash
python3 graph_learning_improvement_experiment.py
```

Requires DeepWalk and reads arc features / edges from CSV.

**4 — Unsupervised U-Net.**
Open the notebook and run top to bottom:

```bash
jupyter notebook Unsup_MSC_UNet/Unsupervised_UNet.ipynb
```

It loads the multi-persistence MSC segmentations as stacked training targets,
trains the U-Net (with a with/without-skip-connections variant), and reports F1
on the test split.

---

## Data

- **DRIVE** — <https://www.isi.uu.nl/Research/Databases/DRIVE/>
- **STARE** — <http://cecas.clemson.edu/~ahoover/stare/>

MSC ground-truth segmentations are expected under
`datasets/.../{drive,stare}_msc_seg/persistence_<p>/` (see the loading cells in
the notebook and the path constants in `mscnn_seg.py`).

---

## Notes & status

This is a 2019 research prototype, published as a study rather than a packaged
tool — expect hard-coded paths, old framework pins, and rough edges. The value
is in the two findings: (1) topology is an excellent *organizer* of an image but
the discriminative signal here was geometric, so neighborhood-similarity graph
learning added little over a random forest on the arc features; and (2)
persistence-varying MSC segmentations are a source of free supervision *and* of
contradiction — managing the latter (via a stability-weighted loss) is the path
to making the unsupervised U-Net work.

## Acknowledgements

The `topoml` package builds on topology tooling by Sam Leventhal, Steve
Petruzza, Attila Gyulassy, Torin McDonald, and Dan Maljovec. GraphSAGE and
DeepWalk are used as vendored components; CAMLP follows Yamaguchi et al.
