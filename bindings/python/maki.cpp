
#include "maki/build/graph/construct_common.hpp"
#include "maki/build/graph/construct_cdbg.hpp"
#include "maki/build/graph/construct_wdbg.hpp"
#include "maki/core/graph/cdbg.hpp"
#include "maki/core/graph/wdbg.hpp"
#include "maki/core/seq/io.hpp"
#include <filesystem>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;

PYBIND11_MODULE(_maki, m) {

  // graph interface

  py::class_<ColouredGraphFiles>(m, "cdbg_files").def(py::init<fs::path>());

  py::class_<ColouredGraph>(m, "cdbg")
      .def(py::init<>())
      .def_static("disk_load", &ColouredGraph::FromDisk,
                  "Load a coloured graph from disk.");

  py::class_<WeightedGraphFiles>(m, "wdbg_files").def(py::init<fs::path>());

  py::class_<WeightedGraph>(m, "wdbg")
      .def(py::init<>())
      .def_static("disk_load", &WeightedGraph::FromDisk,
                  "Load a coloured graph from disk.");

  // build helpers

  py::enum_<InputFileType>(m, "FileType")
      .value("FNA", InputFileType::FastaFileType)
      .value("GFF3", InputFileType::Gff3FileType)
      .value("FQ", InputFileType::FastQFileType)
      .export_values();

  py::class_<GenomeManifest>(m, "genome_manifest")
      .def(py::init<>())
      .def_readonly("files", &GenomeManifest::files)
      .def_readonly("type", &GenomeManifest::type);

  py::class_<DataFilePair>(m, "read_pair")
      .def(py::init<>())
      .def(py::init<std::string&&,std::string&&>(), py::arg("forward"), py::arg("reverse"))
      .def_readwrite("forwardFile", &DataFilePair::forwardFile)
      .def_readwrite("reverseFile", &DataFilePair::reverseFile);

  py::class_<dbg::BuildOptions>(m, "build_opts")
      .def(py::init<>())
      .def(
        py::init<std::size_t,std::size_t,const fs::path&,std::size_t>(),
        py::arg("k"),
        py::arg("s"),
        py::arg("out"),
        py::arg("threads")
      )
      .def_readwrite("kmer_size", &dbg::BuildOptions::kmer_size)
      .def_readwrite("suffix_size", &dbg::BuildOptions::suffix_size)
      .def_readwrite("out", &dbg::BuildOptions::out)
      .def_readwrite("pool_size", &dbg::BuildOptions::pool_size)
      .def_readwrite("reserve_per_chunk", &dbg::BuildOptions::reserve_per_chunk)
      .def_readwrite("threads", &dbg::BuildOptions::threads);

  // build API

  m.def("read_manifest", &readManifest);

  m.def("construct_cdbg",
        py::overload_cast<const GenomeManifest &, dbg::BuildOptions>(
            &cdbg::construct));

  m.def("construct_wdbg",
        py::overload_cast<const DataFilePair &, dbg::BuildOptions>(
            &wdbg::construct));
}
