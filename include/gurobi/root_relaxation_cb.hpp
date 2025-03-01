#ifndef GEDC_ROOT_RELAXATION_CB_HPP
#define GEDC_ROOT_RELAXATION_CB_HPP

#include "gurobi_c++.h"
#include "auxiliary/options.hpp"

class root_relaxation_cb : public GRBCallback {
public:
        options &opt_;
        bool firstvisit_ = true;

        explicit root_relaxation_cb(options &opt) : opt_(opt) {}

protected:
        void callback() override {
                try {
                        if (where != GRB_CB_MIPNODE) {
                                return;
                        }
                        if (getDoubleInfo(GRB_CB_MIPNODE_NODCNT) != 0) {
                                return; // we want root relaxation value
                        }
                        if(firstvisit_) {
                                opt_.ILP_lprelval = getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
                                opt_.ILP_lpreltime = getDoubleInfo(GRB_CB_RUNTIME);
                                firstvisit_ = false;
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

#endif //GEDC_ROOT_RELAXATION_CB_HPP