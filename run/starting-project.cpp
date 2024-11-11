// Copyright © 2021 Giorgio Audrito. All Rights Reserved.


// [INTRODUCTION]
//! Importing the FCPP library.
#include "lib/fcpp.hpp"

//const int BOUND = 25;
const int SOURCE_ID = 0;
const int SINK_ID = 1;

namespace fcpp {

namespace coordination {

namespace tags { 
    struct node_color {};
    struct node_size {};
    struct node_shape {};

    struct is_source {};
    struct is_sink {};

    struct flow_field {};
    struct flow_star_field {};

    //fields of outgoing capacities
    struct capacity_field {};

    struct node_distance_tau{};
    struct node_distance_sigma{};
    struct node_distance_rho{};

    struct out_flow{};
    struct in_flow{};

    struct obstruction_condition{};

    struct dev_num {};
}

//! @brief The maximum communication range between nodes.
constexpr size_t communication_range = 100;

const size_t dim = 2;


FUN field<real_t> capacity_v0(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid!=id ;}, ids);
}


FUN field<real_t> capacity_v1(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid<id ;}, ids);
}



FUN field<real_t> capacity_v2(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid<id ? id-node.uid:node.uid-id;}, ids);
}


FUN field<real_t> capacity_v3(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid<id ? id-node.uid:0;}, ids);
}


FUN field<real_t> capacity(ARGS){ CODE

    return capacity_v2(CALL);
}


real_t sum(field<real_t> input){

    real_t tmp= 0.0;
    for (fcpp::details::field_iterator<field<real_t>> it(input); !it.end(); ++it){
        tmp += it.value();
    }
    return tmp;
}

field<real_t> truncate(field<real_t> inputField, real_t inputValue){

    return inputValue >= 0
    ? map_hood([&](real_t i){
                                        real_t tmp= std::max(std::min(i, inputValue),0.0);
                                        inputValue-= tmp;
                                        return tmp;
                                    }, inputField)
    : -truncate(-inputField, - inputValue);
}

FUN real_t excess(ARGS, field<real_t> flow){ CODE

    bool& is_source_ = node.storage(tags::is_source{});
    bool& is_sink_ = node.storage(tags::is_sink{});


    return is_source_
                ?INF
                :is_sink_
                    ? -INF
                    :sum( flow);
}

FUN int tau_round(ARGS, field<real_t> flow_star){ CODE

    bool& is_sink_ = node.storage(tags::is_sink{});
    bool& is_source_ = node.storage(tags::is_source{});

    field<real_t> graph = capacity(CALL) + flow_star;

    field<bool> is_not_source_field = nbr(CALL, !is_source_);

    return nbr(CALL, 0, [&](field<int> rounds){
            field<real_t> tmp = map_hood([&](int r, real_t g, bool b){
                return g>0 && b ? r : 0;
            }, rounds, graph, is_not_source_field);

            int old_round = self(CALL, rounds);

            int m = max_hood(CALL, tmp);

            return is_sink_
                    ? old_round + 1
                    : std::max(m, old_round);
    });
}




FUN real_t tau(ARGS, field<real_t> flow_star, field<int> tau_round_star){ CODE
    bool& is_sink_ = node.storage(tags::is_sink{});
    bool& is_source_ = node.storage(tags::is_source{});

    field<real_t> graph = capacity(CALL) + flow_star;

    int old_tau_round_ = self(CALL, tau_round_star);

    field<bool> is_not_source_field = nbr(CALL, !is_source_);


    return nbr(CALL, is_sink_? 0.0 : INF, [&](field<real_t> distances){
            field<real_t> tmp = map_hood([&](real_t d, real_t g, int r, bool b){
                return g>0 && r>old_tau_round_ && b 
                        ? d+1 
                        : INF;
            }, distances, graph, tau_round_star, is_not_source_field );

            real_t m = min_hood(CALL, tmp);

            return  is_sink_? 0.0 : m ;
            });
}


FUN real_t sigma(ARGS, field<real_t> flow_star, field<int> tau_round_star){ CODE
    bool& is_source_ = node.storage(tags::is_source{});

    int old_tau_round_ = self(CALL, tau_round_star);

    return nbr(CALL, is_source_? 0.0 : INF, [&](field<real_t> distances){
            field<real_t> tmp = map_hood([&](real_t d, real_t f, real_t r){
                return f>0 && r<= old_tau_round_? d+1 : INF;
            }, distances, flow_star, tau_round_star);

            real_t m = min_hood(CALL, tmp);

            return  is_source_? 0.0 :m;
    });
}

FUN real_t rho(ARGS, field<real_t> flow_star){ CODE
    bool& is_sink_ = node.storage(tags::is_sink{});

    return nbr(CALL, is_sink_? 0.0 : INF, [&](field<real_t> distances){
            field<real_t> tmp = map_hood([&](real_t d, real_t f){
                return f<0 ? d : INF;
            }, distances, flow_star);


            real_t m = min_hood(CALL, tmp) +1;

            return  is_sink_? 0.0 :m;
    });
}



FUN field<real_t> update_flow(ARGS, field<real_t>& flow_star){ CODE
        real_t& tau_= node.storage(tags::node_distance_tau{});
        real_t& sigma_ = node.storage(tags::node_distance_sigma{});
        real_t& rho_ = node.storage(tags::node_distance_rho{});
        field<real_t>& capacity_n = node.storage(tags::capacity_field{});
        field<real_t>& result = node.storage(tags::flow_field{});
        bool& is_source_ = node.storage(tags::is_source{});

        capacity_n = capacity(CALL);

        real_t excess_n = excess(CALL, flow_star);

        int tau_round_ =  tau_round(CALL , flow_star);
        field<int> tau_round_star = nbr(CALL, 0, tau_round_);
        int old_tau_round_ = self(CALL, tau_round_star);

        tau_ = tau(CALL, flow_star, tau_round_star);

        

        sigma_ = sigma(CALL, flow_star, tau_round_star);

        rho_  = rho(CALL, flow_star);

        field<real_t> forward = truncate( (capacity_n + flow_star)
                                            * (nbr(CALL, tau_)< tau_)
                                            * (nbr(CALL, !is_source_)),
                                        excess_n);

        field<real_t> backward = truncate(flow_star 
                                            * (nbr(CALL, sigma_)< sigma_)
                                            * (old_tau_round_ >= tau_round_star),
                                        excess_n);

        field<real_t> reduce = truncate(flow_star
                                            * (nbr(CALL, rho_)< rho_),
                                        excess_n);
        

        result = -flow_star + mux(excess_n<0, 
                                reduce, 
                                mux(tau_<INF, 
                                        forward ,
                                        backward));
        
        return result;
}





//! @brief Main function.
MAIN() {
    // import tag names in the local scope.
    using namespace tags;

    // References
    
    real_t& tau_ = node.storage(node_distance_tau{});
    real_t& out_flow_ = node.storage(out_flow{});
    real_t& in_flow_ = node.storage(in_flow{});
    real_t& obstruction_condition_ = node.storage(obstruction_condition{});
    field<real_t>& capacity_n = node.storage(tags::capacity_field{});
    field<real_t>& flow_star_ = node.storage(tags::flow_star_field{});
    bool& is_source_ = node.storage(is_source{});
    bool& is_sink_ = node.storage(is_sink{});


    if(node.current_time()> 50 && (node.uid % 30 ==2 )){node.terminate();}

    // Usage of node storage
    is_source_ = node.uid== SOURCE_ID;
    is_sink_ = node.uid == SINK_ID;

    node.storage(node_shape{}) = is_source_
                                    ?shape::star 
                                    :is_sink_
                                        ? shape::tetrahedron
                                        : shape::sphere;

    field<real_t> flow_ = nbr(CALL, field<real_t>(0.0),[&](field<real_t> flow_star){
            flow_star_ = flow_star;
            obstruction_condition_ = tau_round(CALL, flow_star);
            return update_flow(CALL, flow_star);
            });
    
    node.storage(node_size{}) = is_sink_ || is_source_ 
                                ? 20
                                : 10;


    out_flow_= is_source_? sum(mux(flow_>0, flow_, 0.0)) : 0.0;
    in_flow_= is_sink_? sum(mux(flow_>0, 0.0, -flow_)) : 0.0;


    //obstruction_condition_ = is_source_? tau_ : 0.0; //sum(flow_)!=-sum(flow_star);


    node.storage(node_color{}) =  is_source_ || (sum(flow_star_)>0 && !is_sink_)
                                    ? color(GREEN)
                                    : sum(flow_star_)<0 || is_sink_
                                        ? color(RED)
                                        : sum(mux(flow_>0, flow_, 0.0))>0
                                            ? color(YELLOW)
                                            : color(BLACK);
}
//! @brief Export types used by the main function.
FUN_EXPORT main_t = export_list<device_t, field<real_t>, real_t, bool, int,  field<int>>;

} // namespace coordination

// [SYSTEM SETUP]

//! @brief Namespace for component options.
namespace option {

//! @brief Import tags to be used for component options.
using namespace component::tags;
//! @brief Import tags used by aggregate functions.
using namespace coordination::tags;

//! @brief Dimensionality of the space.
constexpr size_t dim = 2;

//! @brief Description of the round schedule. To synchronize use:
//using round_s = sequence::periodic_n<1,0,1>;
using round_s = sequence::periodic<
    distribution::interval_n<times_t, 0, 1>,    // uniform time in the [0,1] interval for start
    distribution::weibull_n<times_t, 10, 1, 10> // weibull-distributed time for interval (10/10=1 mean, 1/10=0.1 deviation)
>;

//! @brief The sequence of network snapshots (one every simulated second).
using log_s = sequence::periodic_n<1, 0, 1>;
//! @brief The sequence of node generation events (dev_num devices all generated at time 0).
using spawn_s = sequence::multiple<distribution::constant_i<size_t,dev_num>, distribution::constant_n<real_t,0>>;
//! @brief The distribution of initial node positions.
using rectangle_d = distribution::rect_n<1, 0, 0, 1000, 1000>;
//! @brief The contents of the node storage as tags and associated types. 
using store_t = tuple_store<
    node_color,                         color,
    node_size,                          double,
    node_shape,                         shape,
    node_distance_tau,                  real_t,
    node_distance_sigma,                real_t,
    node_distance_rho,                   real_t,
    capacity_field,                     field<real_t>,
    flow_field,                         field<real_t>,
    flow_star_field,                    field<real_t>,
    out_flow,                           real_t,
    in_flow,                            real_t,
    obstruction_condition,              real_t,
    is_source,                          bool,
    is_sink,                            bool
>;
//! @brief The tags and corresponding aggregators to be logged (change as needed).
using aggregator_t = aggregators< 
    out_flow,                   aggregator::max<real_t>,
    in_flow,                    aggregator::max<real_t>,
    obstruction_condition,      aggregator::min<real_t>
>;
using plot_t = plot::split<
    dev_num,
    plot::split<
        seed,
        plot::split<
            plot::time,
            plot::join<
                plot::value<aggregator::max<in_flow>>,
                plot::value<aggregator::max<out_flow>>
            >
        >
    >
>;

//! @brief The general simulation options.
DECLARE_OPTIONS(list,
    parallel<true>,      // multithreading enabled on node rounds
    synchronised<false>, // optimise for asynchronous networks
    program<coordination::main>,   // program to be run (refers to MAIN above)
    exports<coordination::main_t>, // export type list (types used in messages)
    retain<metric::retain<2,1>>,   // messages are kept for 2 seconds before expiring
    round_schedule<round_s>, // the sequence generator for round events on nodes
    log_schedule<log_s>,     // the sequence generator for log events on the network
    spawn_schedule<spawn_s>, // the sequence generator of node creation events on the network
    store_t,       // the contents of the node storage
    extra_info<dev_num, int, seed, int>,
    aggregator_t,  // the tags and corresponding aggregators to be logged
    plot_type<plot_t>,
    init<
        x,      rectangle_d // initialise position randomly in a rectangle for new nodes
    >,
    dimension<dim>, // dimensionality of the space
    connector<connect::fixed<100, 1, dim>>, // connection allowed within a fixed comm range
    //label_text_tag<node_text>, // the text to associate to a node
    label_size_val<1, 4>,        // the size of the text
    shape_tag<node_shape>, // the shape of a node is read from this tag in the store
    size_tag<node_size>,   // the size  of a node is read from this tag in the store
    color_tag<node_color>  // the color of a node is read from this tag in the store
);

} // namespace option

} // namespace fcpp


//! @brief The main function.
int main() {
    using namespace fcpp;

    // The plotter object.
    option::plot_t plotter;
    // The network object type (interactive simulator with given options).
    using net_t = component::interactive_simulator<option::list>::net;
    std::cout << "/*\n";
    for (int seed=6; seed<7; seed++) for (int num=600; num<=900; num+=300) {
        // The initialisation values (simulation name).
        auto init_v = common::make_tagged_tuple<option::name, option::dev_num, option::seed, option::plotter>("Starting Project", num, seed, &plotter);
        // Construct the network object.
        net_t network{init_v};
        // Run the simulation until exit.
        network.run();
    }
    // Builds the resulting plots.
    std::cout << "*/\n" << plot::file("starting-project", plotter.build());
    return 0;
}
