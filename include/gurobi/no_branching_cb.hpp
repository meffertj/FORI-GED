#ifndef GEDC_NO_BRANCHING_CB_HPP
#define GEDC_NO_BRANCHING_CB_HPP

#include "gurobi_c++.h"
#include "auxiliary/options.hpp"

class no_branching_cb : public GRBCallback {
public:
        options &opt_;

        explicit no_branching_cb(options &opt) : opt_(opt) {}

protected:
        void callback() override {
                try {
                        if (where != GRB_CB_MIPNODE) {
                                return;
                        }
                        if (getDoubleInfo(GRB_CB_MIPNODE_NODCNT) != 0) {
                                abort();
                        }
                        opt_.root_primal_ = getDoubleInfo(GRB_CB_MIPNODE_OBJBST);
                        opt_.root_dual_ = getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
                        opt_.root_time_ = getDoubleInfo(GRB_CB_RUNTIME);
                }
                catch (GRBException &e) {
                        std::cout << "Error number: " << e.getErrorCode() << std::endl;
                        std::cout << e.getMessage() << std::endl;
                }
                catch (...) {
                        std::cout << "Error during callback" << std::endl;
                }
        }
};

#endif //GEDC_NO_BRANCHING_CB_HPP
