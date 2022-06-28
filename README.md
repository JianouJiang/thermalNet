<div align="center">
  <img src="https://github.com/JianouJiang/thermalNet/blob/main/thermalNet_logo.png">
</div>

[![Python](https://img.shields.io/pypi/pyversions/tensorflow.svg?style=plastic)](https://badge.fury.io/py/tensorflow)
[![PyPI](https://badge.fury.io/py/tensorflow.svg)](https://badge.fury.io/py/tensorflow)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4724125.svg)](https://doi.org/10.5281/zenodo.4724125)

**`Documentation`** |
------------------- |
[![Documentation](https://img.shields.io/badge/api-reference-blue.svg)](https://github.com/JianouJiang/thermalNet/README.md) |

[thermalNet](https://www.thermalnet.org/) is an end-to-end open source physics-informed
 neural network platform based on [TensorFlow](https://www.tensorflow.org/) for solving 
 heat transfer related problems. It has a comprehensive, flexible ecosystem of
[examples](https://github.com/JianouJiang/thermalNet/examples),
[tools](https://github.com/JianouJiang/thermalNet/tools), and
[community](https://www.thermalnet.org/community) resources that lets
researchers and developers pull the framework of this PINN to validate and 
upgrade their optimised numerical methods.

ThermalNet was originally developed by researchers and engineers working on the
[osneyNet](https://github.com/osneyNet) team within top UK universities to conduct a physically constrained
 deep neural networks research. The system is general enough to be applicable 
 in a other physics-informed problems beyond heat transfer problem too.

ThermalNet provides stable [Python](https://www.thermalnet.org/api_docs/python)
APIs, as well as non-guaranteed backward compatible API for
[other languages](https://www.thermalnet.org/api_docs).


## Install

See the [thermalNet install guide](https://www.thermalnet.org/install) for the
[pip package](https://www.thermalnet.org/install/pip), to
[enable GPU support](https://www.thermalnet.org/install/gpu), use a
[Docker container](https://www.thermalnet.org/install/docker), and
[build from source](https://www.thermalnet.org/install/source).

To install the current release, which includes support for
[CUDA-enabled GPU cards](https://www.thermalnet.org/install/gpu) *(Ubuntu and
Windows)*:

```
$ pip install thermalnet
```

To update thermalNet to the latest version, add `--upgrade` flag to the above
commands.

*Nightly binaries are available for testing using the
[tf-nightly](https://pypi.python.org/pypi/tf-nightly) and
[tf-nightly-cpu](https://pypi.python.org/pypi/tf-nightly-cpu) packages on PyPi.*

#### *Try your first thermalNet program*

```shell
$ python
```

```python
>>> import thermalnet as tn
>>> a_list = [1,9,4,9,1,0,0,1]
>>> tn.normalisation(a_list)
[0.11, 1.0, 0.44, 1.0, 0.11, 0.0, 0.0, 1.0]
>>> tn.test()
'Contributors introduction...'
```

For more examples, see the
[thermalNet tutorials](https://www.thermalnet.org/tutorials/).

## Contribution guidelines

**If you want to contribute to thermalNet, be sure to review the
[contribution guidelines](CONTRIBUTING.md). This project adheres to thermalNet's
[code of conduct](CODE_OF_CONDUCT.md). By participating, you are expected to
uphold this code.**

**We use [GitHub issues](https://github.com/thermalnet/thermalnet/issues) for
tracking requests and bugs, please see
[thermalNet Discuss](https://groups.google.com/a/thermalnet.org/forum/#!forum/discuss)
for general questions and discussion, and please direct specific questions to
[Stack Overflow](https://stackoverflow.com/questions/tagged/thermalnet).**

The ThermalNet project strives to abide by generally accepted best practices in
open-source software development:

[![Fuzzing Status](https://oss-fuzz-build-logs.storage.googleapis.com/badges/tensorflow.svg)](https://bugs.chromium.org/p/oss-fuzz/issues/list?sort=-opened&can=1&q=proj:tensorflow)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/1486/badge)](https://bestpractices.coreinfrastructure.org/projects/1486)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v1.4%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)



## Citation

    @article{jianou2022physics,
      title={Physics-informed neural networks: n/a},
      journal={Journal of Computational Physics},
      volume={n/a},
      pages={n/a-n/a},
      year={2022},
      publisher={Elsevier}
    }

  - Jianou et al. "[Physics-informed neural networks: n/a](https://www.sciencedirect.com/science/article/pii/n/a)." Journal of Computational Physics n/a (2022): n/a-n/a.

### Official Builds

Build Type                    | Status                                                                                                                                                                           | Artifacts
----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------
**Linux CPU**                 | [![Status](https://storage.googleapis.com/tensorflow-kokoro-build-badges/ubuntu-cc.svg)](https://storage.googleapis.com/tensorflow-kokoro-build-badges/ubuntu-cc.html)           | [PyPI](https://pypi.org/project/tf-nightly/)
**Linux GPU**                 | [![Status](https://storage.googleapis.com/tensorflow-kokoro-build-badges/ubuntu-gpu-py3.svg)](https://storage.googleapis.com/tensorflow-kokoro-build-badges/ubuntu-gpu-py3.html) | [PyPI](https://pypi.org/project/tf-nightly-gpu/)
**Linux XLA**                 | [![Status](https://storage.googleapis.com/tensorflow-kokoro-build-badges/ubuntu-xla.svg)](https://storage.googleapis.com/tensorflow-kokoro-build-badges/ubuntu-xla.html)         | TBA
**macOS**                     | [![Status](https://storage.googleapis.com/tensorflow-kokoro-build-badges/macos-py2-cc.svg)](https://storage.googleapis.com/tensorflow-kokoro-build-badges/macos-py2-cc.html)     | [PyPI](https://pypi.org/project/tf-nightly/)
**Windows CPU**               | [![Status](https://storage.googleapis.com/tensorflow-kokoro-build-badges/windows-cpu.svg)](https://storage.googleapis.com/tensorflow-kokoro-build-badges/windows-cpu.html)       | [PyPI](https://pypi.org/project/tf-nightly/)
**Windows GPU**               | [![Status](https://storage.googleapis.com/tensorflow-kokoro-build-badges/windows-gpu.svg)](https://storage.googleapis.com/tensorflow-kokoro-build-badges/windows-gpu.html)       | [PyPI](https://pypi.org/project/tf-nightly-gpu/)
**Android**                   | [![Status](https://storage.googleapis.com/tensorflow-kokoro-build-badges/android.svg)](https://storage.googleapis.com/tensorflow-kokoro-build-badges/android.html)               | [Download](https://bintray.com/google/tensorflow/tensorflow/_latestVersion)
**Raspberry Pi 0 and 1**      | [![Status](https://storage.googleapis.com/tensorflow-kokoro-build-badges/rpi01-py3.svg)](https://storage.googleapis.com/tensorflow-kokoro-build-badges/rpi01-py3.html)           | [Py3](https://storage.googleapis.com/tensorflow-nightly/tensorflow-1.10.0-cp34-none-linux_armv6l.whl)
**Raspberry Pi 2 and 3**      | [![Status](https://storage.googleapis.com/tensorflow-kokoro-build-badges/rpi23-py3.svg)](https://storage.googleapis.com/tensorflow-kokoro-build-badges/rpi23-py3.html)           | [Py3](https://storage.googleapis.com/tensorflow-nightly/tensorflow-1.10.0-cp34-none-linux_armv7l.whl)
**Libtensorflow MacOS CPU**   | Status Temporarily Unavailable                                                                                                                                                   | [Nightly Binary](https://storage.googleapis.com/libtensorflow-nightly/prod/tensorflow/release/macos/latest/macos_cpu_libtensorflow_binaries.tar.gz) [Official GCS](https://storage.googleapis.com/tensorflow/)
**Libtensorflow Linux CPU**   | Status Temporarily Unavailable                                                                                                                                                   | [Nightly Binary](https://storage.googleapis.com/libtensorflow-nightly/prod/tensorflow/release/ubuntu_16/latest/cpu/ubuntu_cpu_libtensorflow_binaries.tar.gz) [Official GCS](https://storage.googleapis.com/tensorflow/)
**Libtensorflow Linux GPU**   | Status Temporarily Unavailable                                                                                                                                                   | [Nightly Binary](https://storage.googleapis.com/libtensorflow-nightly/prod/tensorflow/release/ubuntu_16/latest/gpu/ubuntu_gpu_libtensorflow_binaries.tar.gz) [Official GCS](https://storage.googleapis.com/tensorflow/)
**Libtensorflow Windows CPU** | Status Temporarily Unavailable                                                                                                                                                   | [Nightly Binary](https://storage.googleapis.com/libtensorflow-nightly/prod/tensorflow/release/windows/latest/cpu/windows_cpu_libtensorflow_binaries.tar.gz) [Official GCS](https://storage.googleapis.com/tensorflow/)
**Libtensorflow Windows GPU** | Status Temporarily Unavailable                                                                                                                                                   | [Nightly Binary](https://storage.googleapis.com/libtensorflow-nightly/prod/tensorflow/release/windows/latest/gpu/windows_gpu_libtensorflow_binaries.tar.gz) [Official GCS](https://storage.googleapis.com/tensorflow/)

## Resources

*   [thermalNet.org](https://www.thermalnet.org)
*   [thermalNet Tutorials](https://www.thermalNet.org/tutorials/)
*   [thermalNet Examples](https://github.com/thermalnet/examples)
*   [DeepLearning.AI TensorFlow Developer Professional Certificate](https://www.coursera.org/specializations/tensorflow-in-practice)
*   [TensorFlow: Data and Deployment from Coursera](https://www.coursera.org/specializations/tensorflow-data-and-deployment)
*   [Getting Started with TensorFlow 2 from Coursera](https://www.coursera.org/learn/getting-started-with-tensor-flow2)
*   [TensorFlow: Advanced Techniques from Coursera](https://www.coursera.org/specializations/tensorflow-advanced-techniques)
*   [TensorFlow 0.12.0](https://storage.googleapis.com/tensorflow/windows/cpu/tensorflow-0.12.0rc0-cp35-cp35m-win_amd64.whl)
*   [Physics Informed Neural Networks](https://github.com/maziarraissi/PINNs)
*   [Python 3.5.4](https://www.python.org/downloads/release/python-354/)
*   [2D heat diffusion equation using Python](https://www.youtube.com/watch?v=mSYm46VVZRo&t=17s)
*   [BVPs for Heat Equations (Exact Solutions)](https://www.cfm.brown.edu/people/dobrush/am34/Mathematica/ch6/bheat.html)
*   [Finite difference schemes for multilayer diffusion](https://www.sciencedirect.com/science/article/pii/S0895717711000938)
*   [Multigrid solver for 2D heat conduction problems](https://www.researchgate.net/publication/334774068_Multigrid_solver_for_2D_heat_conduction_problems)

Learn more about the
[ThermalNet community](https://www.thermalnet.org/community) and how to
[contribute](https://www.thermalnet.org/community/contribute).

## License

[Apache License 2.0](LICENSE)
