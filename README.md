# SDDPBlock

SMS++ Block and Solver for multistage stochastic programming problems.

## Getting started

These instructions will let you build SDDPBlock on your system.

### Requirements

- [SMS++ StochasticBlock](https://gitlab.com/smspp/stochasticblock)
- [STochastic OPTimization library (StOpt)](https://gitlab.com/stochastic-control/StOpt)

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
make
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
sudo make install
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(SDDPBlock)
target_link_libraries(<my_target> SMS++::SDDPBlock)
```

## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/sddpblock/-/issues/new).

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.

## Authors

### Current Lead Authors

- **Rafael Durbano Lobato**  
  Dipartimento di Informatica  
  Università di Pisa

## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.

## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
