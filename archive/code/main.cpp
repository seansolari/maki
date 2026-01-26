#include <map>
#include <string>
#include <string_view>
#include <glog/logging.h>
#include <maki/maki.h>
#include <maki/cli.hpp>
#include <maki/version.hpp>

int main(int argc, char const *argv[]) {
    seqan3::argument_parser base_parser{
        "maki",
        argc,
        argv,
        seqan3::update_notifications::off,
        {"index", "dbstats", "make-filter", "classify", "dist", "count", "subcount", "jaccard",
         "compress1", "lca1", "glca1", "randist", "grep", "coverage", "kmers"}};

    // top-level metadata
    base_parser.info.author = "sean.solari@monash.edu";
    base_parser.info.description.push_back("Perform indexing, distance calculations or classification from genome sequences.");
    base_parser.info.version = MAKI::MAKI_version_cstring;

    // parse command
    try { base_parser.parse(); }
    catch (seqan3::argument_parser_error const &ext) { LOG(ERROR) << ext.what(); return -1; }
    seqan3::argument_parser &command_parser = base_parser.get_sub_parser();
    
    std::unique_ptr<BaseCommand> mainfp;
         if (command_parser.info.app_name == "maki-index")          mainfp = std::make_unique<IndexCommand>();
    else if (command_parser.info.app_name == "maki-dbstats")        mainfp = std::make_unique<DbStatsCommand>();
    else if (command_parser.info.app_name == "maki-make-filter")    mainfp = std::make_unique<MakeFilterCommand>();
    else if (command_parser.info.app_name == "maki-classify")       mainfp = std::make_unique<ClassifyCommand>();
    else if (command_parser.info.app_name == "maki-dist")           mainfp = std::make_unique<DistanceCommand>();
    else if (command_parser.info.app_name == "maki-count")          mainfp = std::make_unique<CountCommand>();
    else if (command_parser.info.app_name == "maki-coverage")       mainfp = std::make_unique<CoverageCommand>();
    else if (command_parser.info.app_name == "maki-subcount")       mainfp = std::make_unique<SubcountCommand>();
    else if (command_parser.info.app_name == "maki-jaccard")        mainfp = std::make_unique<JaccardCommand>();
    else if (command_parser.info.app_name == "maki-compress1")      mainfp = std::make_unique<CompressCommand>();
    else if (command_parser.info.app_name == "maki-lca1")           mainfp = std::make_unique<LcaCommand>();
    else if (command_parser.info.app_name == "maki-glca1")          mainfp = std::make_unique<GlcaCommand>();
    else if (command_parser.info.app_name == "maki-randist")        mainfp = std::make_unique<SamdistCommand>();
    else if (command_parser.info.app_name == "maki-grep")           mainfp = std::make_unique<QueryCommand>();
    else if (command_parser.info.app_name == "maki-kmers")           mainfp = std::make_unique<KmersCommand>();
    else {
        LOG(ERROR) << "unrecognised command " << command_parser.info.app_name;
        return -1;
    }
    return mainfp->run(argv[0], command_parser);
}
