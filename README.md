---
languages:
- Matlab
description: "Demo code for encoding perceptually salient early reflections for parametric spatial audio rendering"
urlFragment: "https://github.com/microsoft/Perceptual_saliency_of_early_reflections"
---

# Official Microsoft Sample

Demo code for encoding perceptually salient early reflections for parametric spatial audio rendering.

## Contents

Outline the file contents of the repository. It helps users navigate the codebase, build configuration and any related assets.

| File/folder       | Description                                |
|-------------------|--------------------------------------------|
| `demo.m`          | Run this file to start the demo.           |
| `Code/`           | Functions to run `demo.m`                  |
| `Data/`           | Data required by `demo.m`                  |
|                   |                                            |
| `playAudio.pd`    | Use this file to listen to audio examples, |
|                   | or manually listen to the wav files.       |
| `Auralizations/`  | Audio examples used by `playAudio.pd`      |
|                   |                                            |
| `.gitignore`      | Define what to ignore at commit time.      |
| `README.md`       | This README file.                          |
| `LICENSE`         | The license for the sample.                |

## Prerequisites

Needs parts of the FABIAN HRTF data base from (https://dx.doi.org/10.14279/depositonce-5718.3). Data can be loaded from `demo.m`.

## Setup

Run `demo.m` (will tmporarily add needed folders to the Matlab search path).

## Running the sample

Run `demo.m` or execute block by block.
Listen to audio files in `Auralizations/` or use `playAudio.pd`.

## Contributing

This project welcomes contributions and suggestions.  Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit https://cla.opensource.microsoft.com.

When you submit a pull request, a CLA bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., status check, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.
