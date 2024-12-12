import contagion
import create_network


def main():

    config = contagion.SimpleContagionSim.checked_config_load("config.yaml")

    network = create_network.create_twolayer_bipartite_network(
        n_indvs=100,
        n_hh=40,
        n_comm=2,
        p_comm=0.4,
    )

    initial_condition = contagion.SimpleContagionSim.create_initial_state(
        graph=network,
        n_seeds=config["seeding"]["num"],
    )

    sim = contagion.SimpleContagionSim(
        graph=network,
        initial_state=initial_condition,
        parameters=config["parameters"],
        return_statuses="SEIRDTQ",
    )

    sim.run()
    sim.write("sample_results.h5")


if __name__ == "__main__":
    main()
