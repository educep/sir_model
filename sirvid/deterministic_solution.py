"""
Created by ecepeda at 19/11/2021
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use("bmh")
color_ = sns.color_palette("Paired")


def sirvid_model(
    days,
    kappa,
    n_population,
    s_init_perc,
    vac_init_perc,
    recovery_rate,
    infection_rate,
    vaccination_rate,
    recovery_rate_vac,
    vac_protect,
    inf_init_perc,
    inf_vac_init_perc,
):
    """
    defining parameters
    """
    # days = 15  # average duration of disease
    # kappa = 2.5  # avg nb of interactions by day by individual
    # n_population = 1e6  # population size
    # s_init_perc = .985  # S0 (initial susceptible set as % of pop)
    # vac_init_perc = .005  # V0 (initial vaccinated set as % of pop)
    # inf_vac_init_perc = vac_init_perc * 0.1
    # inf_init_perc = 1 - s_init_perc - vac_init_perc - inf_vac_init_perc  # I0 (initial infected set as % of pop)

    # parameters for non vaccinated people
    # recovery_rate = .90  # recovery rate from disease
    mortality = 1 - recovery_rate  # rate decease from disease
    # infection_rate = 0.14  # probability of infection for one interaction between susceptible and infected
    beta = kappa * infection_rate  # transmission rate
    gamma = recovery_rate / days  # recovery rate by day
    mu = mortality / days  # decease rate by day

    # parameters for vaccinated people
    # vaccination_rate = 0.005  # vaccination rate
    # recovery_rate_vac = .99  # recovery rate from disease
    mortality_vac = 1 - recovery_rate_vac  # rate decease from disease
    # vac_protect = .95  # decrease of infection rate thanks to vaccine as %
    infection_rate_vac = infection_rate * (
        1 - vac_protect
    )  # infection proba between vaccinated and infected
    beta_vac = kappa * infection_rate_vac  # transmission rate for vaccinated people
    gamma_vac = recovery_rate_vac / days  # recovery rate by day
    mu_vac = mortality_vac / days  # decease rate by day
    nbdays = 200  # total number of days for simulation

    print(
        "\ninitial susceptible pct = {:.2f}%\ninitial vaccinated pct = {:.2f}%\n"
        "initial susceptible infected pct = {:.2f}%\ninitial vaccinated infected pct = {:.2f}%".format(
            s_init_perc * 100,
            vac_init_perc * 100,
            inf_init_perc * 100,
            inf_vac_init_perc * 100,
        )
    )

    assert (
        abs(inf_init_perc + inf_vac_init_perc + s_init_perc + vac_init_perc - 1) < 1e-8
    )

    aa = (
        "\nrecup rate = {:.2f}% vaccinated = {:.2f}%\n"
        "mortality rate = {:.2f}% vaccinated = {:.2f}%\n"
        "proba infection = {:.2f}% vaccinated = {:.2f}%\n"
        "beta = {:.2f} vaccinated = {:.2f}\n"
        "gamma = {:.2f} vaccinated = {:.2f}"
    )

    print(
        aa.format(
            recovery_rate * 100,
            recovery_rate_vac * 100,
            mu * 100,
            mu_vac * 100,
            infection_rate * 100,
            infection_rate_vac * 100,
            beta,
            beta_vac,
            gamma,
            gamma_vac,
        )
    )

    # declaring objects
    time_frame = range(0, nbdays)
    data = pd.DataFrame(
        0,
        index=time_frame,
        columns=[
            "susceptibles",
            "infected",
            "recovered",
            "deceased",
            "vaccinated",
            "infected_vac",
            "recovered_vac",
            "deceased_vac",
        ],
    )
    increments = pd.DataFrame(
        0,
        index=time_frame,
        columns=[
            "vaccinations",
            "infections",
            "recoveries",
            "deceases",
            "infections_vac",
            "recoveries_vac",
            "deceases_vac",
        ],
    )

    # initializing instances, the other population compartments are set to 0 (always)
    # problem, if vaccinated people are significant, say > 10%, we should consider a % > 0 of infected_vac
    data.loc[0, "susceptibles"] = n_population * s_init_perc
    data.loc[0, "vaccinated"] = n_population * vac_init_perc
    data.loc[0, "infected"] = n_population * inf_init_perc
    data.loc[0, "infected_vac"] = n_population * inf_vac_init_perc

    for day_ in range(1, nbdays):
        # day_ = 1
        prev_day = day_ - 1
        # computing increments for non vaccinated people
        suscept_prev_day = data.loc[prev_day, "susceptibles"]
        infected_prev_day = data.loc[prev_day, "infected"]
        total_infected_prev_day = infected_prev_day + data.loc[prev_day, "infected_vac"]
        alive_pop_prev_day = data.loc[
            prev_day,
            [
                "susceptibles",
                "infected",
                "recovered",
                "vaccinated",
                "infected_vac",
                "recovered_vac",
            ],
        ].sum()

        increments.loc[day_, "vaccinations"] = vaccination_rate * suscept_prev_day
        increments.loc[day_, "infections"] = (
            beta * suscept_prev_day * total_infected_prev_day / alive_pop_prev_day
        )
        increments.loc[day_, "recoveries"] = gamma * infected_prev_day
        increments.loc[day_, "deceases"] = mu * infected_prev_day

        # computing increments for vaccinated people
        vac_prev_day = data.loc[prev_day, "vaccinated"]
        infected_vac_prev_day = data.loc[prev_day, "infected_vac"]
        increments.loc[day_, "infections_vac"] = (
            beta_vac * vac_prev_day * total_infected_prev_day / alive_pop_prev_day
        )
        increments.loc[day_, "recoveries_vac"] = gamma_vac * infected_vac_prev_day
        increments.loc[day_, "deceases_vac"] = mu_vac * infected_vac_prev_day

        # aggregating increments for non vaccinated people
        data.loc[day_, "susceptibles"] = (
            suscept_prev_day
            - increments.loc[day_, ["vaccinations", "infections"]].sum()
        )
        data.loc[day_, "infected"] = (
            infected_prev_day
            + increments.loc[day_, "infections"]
            - increments.loc[day_, ["recoveries", "deceases"]].sum()
        )
        data.loc[day_, "recovered"] = (
            data.loc[prev_day, "recovered"] + increments.loc[day_, "recoveries"]
        )
        data.loc[day_, "deceased"] = (
            data.loc[prev_day, "deceased"] + increments.loc[day_, "deceases"]
        )

        # aggregating increments for vaccinated people
        data.loc[day_, "vaccinated"] = (
            data.loc[prev_day, "vaccinated"]
            + increments.loc[day_, "vaccinations"]
            - increments.loc[day_, "infections_vac"]
        )
        data.loc[day_, "infected_vac"] = (
            data.loc[prev_day, "infected_vac"]
            + increments.loc[day_, "infections_vac"]
            - increments.loc[day_, ["recoveries_vac", "deceases_vac"]].sum()
        )
        data.loc[day_, "recovered_vac"] = (
            data.loc[prev_day, "recovered_vac"] + increments.loc[day_, "recoveries_vac"]
        )
        data.loc[day_, "deceased_vac"] = (
            data.loc[prev_day, "deceased_vac"] + increments.loc[day_, "deceases_vac"]
        )

    return data, increments


def do_summary(data, inc):
    n_population_ = int(data.iloc[0].sum())
    summary = pd.DataFrame(columns=["No Vacunados", "Vacunados", "Total"])
    summary.loc["Pob Vacunada", "No Vacunados"] = 0
    summary.loc["Pob Vacunada", "Vacunados"] = (
        inc["vaccinations"].sum() + data["vaccinated"].iloc[0]
    )
    summary.loc["Contagios", "No Vacunados"] = (
        inc["infections"].sum() + data["infected"].iloc[0]
    )
    summary.loc["Contagios", "Vacunados"] = (
        inc["infections_vac"].sum() + data["infected_vac"].iloc[0]
    )
    summary.loc["Recuperaciones", "No Vacunados"] = (
        inc["recoveries"].sum() + data["recovered"].iloc[0]
    )
    summary.loc["Recuperaciones", "Vacunados"] = (
        inc["recoveries_vac"].sum() + data["recovered_vac"].iloc[0]
    )
    summary.loc["Decesos", "No Vacunados"] = (
        inc["deceases"].sum() + data["deceased"].iloc[0]
    )
    summary.loc["Decesos", "Vacunados"] = (
        inc["deceases_vac"].sum() + data["deceased_vac"].iloc[0]
    )
    summary.loc["No Afectados", "Vacunados"] = (
        summary.loc["Pob Vacunada", "Vacunados"] - summary.loc["Contagios", "Vacunados"]
    )
    total_unaffected = n_population_ - summary.loc["Contagios"].sum()
    summary.loc["No Afectados", "No Vacunados"] = (
        total_unaffected - summary.loc["No Afectados", "Vacunados"]
    )

    summary.loc["Max Nb Infectados", "No Vacunados"] = data["infected"].max()
    summary.loc["Max Nb Infectados", "Vacunados"] = data["infected_vac"].max()

    summary.loc["Max Nb Contagios diarios", "No Vacunados"] = inc["infections"].max()
    summary.loc["Max Nb Contagios diarios", "Vacunados"] = inc["infections_vac"].max()
    summary.loc["Max Nb Recuperaciones diarias", "No Vacunados"] = inc[
        "recoveries"
    ].max()
    summary.loc["Max Nb Recuperaciones diarias", "Vacunados"] = inc[
        "recoveries_vac"
    ].max()
    summary.loc["Max Nb Decesos diarios", "No Vacunados"] = inc["deceases"].max()
    summary.loc["Max Nb Decesos diarios", "Vacunados"] = inc["deceases_vac"].max()
    summary["Total"] = summary[["No Vacunados", "Vacunados"]].sum(axis=1)

    return summary.astype(int)


if __name__ == "__main__":
    """
    defining parameters SIRVID
    """
    days = 15  # average duration of disease
    kappa = 5  # avg nb of interactions by day by individual
    n_population = 1e6  # population size
    vac_init_perc = 0.005  # V0 (initial vaccinated set as % of pop)
    inf_init_perc = 0.0075  # I0 (initial infected set as % of pop)
    inf_vac_init_perc = 0.0025  # I0 (initial vaccinated and infected set as % of pop)
    s_init_perc = 1 - vac_init_perc - inf_init_perc - inf_vac_init_perc
    #
    # parameters for non vaccinated people
    recovery_rate = 0.90  # recovery rate from disease
    # mortality = 1 - recovery_rate  # rate decease from disease
    infection_rate = 0.1  # probability of infection for one interaction between susceptible and infected
    # beta = kappa * infection_rate  # transmission rate
    # gamma = recovery_rate / days  # recovery rate by day
    # mu = mortality / days  # decease rate by day
    #
    # parameters for vaccinated people
    vaccination_rate = 0.03  # vaccination rate
    recovery_rate_vac = 0.99  # recovery rate from disease
    # mortality_vac = 1 - recovery_rate_vac  # rate decease from disease
    vac_protect = 0.90  # decrease of infection rate thanks to the vaccine as %
    # infection_rate_vac = infection_rate * (1 - vac_protect)  # infection proba between vaccinated and infected
    # beta_vac = kappa * infection_rate_vac  # transmission rate for vaccinated people
    # gamma_vac = recovery_rate_vac / days  # recovery rate by day
    # mu_vac = mortality_vac / days  # decease rate by day
    # nbdays = 150  # total number of days for simulation

    data, inc = sirvid_model(
        days,
        kappa,
        n_population,
        s_init_perc,
        vac_init_perc,
        recovery_rate,
        infection_rate,
        vaccination_rate,
        recovery_rate_vac,
        vac_protect,
        inf_init_perc,
        inf_vac_init_perc,
    )

    summary = do_summary(data, inc)
    summary.to_clipboard()
    # charts
    cols = [
        "susceptibles",
        "vaccinated",
        "infected",
        "infected_vac",
        "recovered",
        "recovered_vac",
        "deceased",
        "deceased_vac",
    ]

    renamed = {
        "vaccinated": "vacunados",
        "infected": "infectados",
        "infected_vac": "vacunados e infectados",
        "recovered": "recuperados",
        "recovered_vac": "vacunados y recuperados",
        "deceased": "decesos",
        "deceased_vac": "vacunas decesos",
    }

    data_ = data[cols[::-1]].rename(columns=renamed)
    data_norm = data_ / n_population

    max_x = 120
    f, ax = plt.subplots(1, 2, figsize=(17, 5))
    ax[0].set_xlim([0, max_x])
    ax[0].set_ylim([0, n_population])
    data_.plot(legend=False, ax=ax[0], color=color_)

    ax[1].set_xlim([0, max_x])
    ax[1].set_ylim([0, 1])
    data_norm.plot.area(legend=False, ax=ax[1], color=color_)

    # box = ax[0].get_position()
    # ax[0].set_position([box.x0 - box.width * 0.02, box.y0, box.width * 0.98, box.height])
    # ax[0].legend(loc='upper center', bbox_to_anchor=(1.2, 0.8), fancybox=False, shadow=False, ncol=1)

    box = ax[1].get_position()
    ax[1].set_position([box.x0 - box.width * 0.1, box.y0, box.width * 0.98, box.height])
    ax[1].legend(
        loc="upper center",
        bbox_to_anchor=(1.2, 0.8),
        fancybox=False,
        shadow=False,
        ncol=1,
    )

    infection_rate = 0.15
    vac_protect = 0.9
    vac_init_perc = [0.005, 0.5, 0.6, 0.8]
    results = {}
    for vac_p in vac_init_perc:
        s_init_perc = 1 - vac_p - inf_init_perc - inf_vac_init_perc
        results[int(100 * vac_p)] = sirvid_model(
            days,
            kappa,
            n_population,
            s_init_perc,
            vac_p,
            recovery_rate,
            infection_rate,
            vaccination_rate,
            recovery_rate_vac,
            vac_protect,
            inf_init_perc,
            inf_vac_init_perc,
        )

    infections = pd.concat(
        [aa[0][["infected", "infected_vac"]].sum(axis=1) for ii, aa in results.items()],
        axis=1,
    )
    infections.columns = ["V0 {}%".format(_) for _ in results.keys()]
    infections = 100 * infections / n_population
    ax = infections.iloc[:150].plot(figsize=(6, 4))
    ax.set_ylim([0, 60])
    ax.set_xlim([0, 150])

    deceases = pd.concat(
        [aa[0][["deceased", "deceased_vac"]].sum(axis=1) for ii, aa in results.items()],
        axis=1,
    )
    deceases.columns = ["V0 {}%".format(_) for _ in results.keys()]
    deceases = 100 * deceases / n_population
    ax = deceases.iloc[:150].plot(figsize=(6, 4), legend=False)
    ax.set_ylim([0, 10])
    ax.set_xlim([0, 150])

    f, ax = plt.subplots(1, 4, figsize=(14, 3))
    for kk, (ii, aa) in enumerate(results.items()):
        deceases_ = 100 * aa[0][["infected", "infected_vac"]] / n_population
        deceases_.columns = ["contagios", "contagios vacunados"]
        deceases_.plot.area(
            ax=ax[kk], legend=kk == 3, color=[color_[i] for i in [5, 4]]
        )
        ax[kk].set_ylim([0, 60])
        ax[kk].set_xlim([0, 150])

    f, ax = plt.subplots(1, 4, figsize=(14, 3))
    for kk, (ii, aa) in enumerate(results.items()):
        deceases_ = 100 * aa[0][["deceased", "deceased_vac"]] / n_population
        deceases_.columns = ["decesos", "decesos vacunados"]
        deceases_.plot.area(
            ax=ax[kk], legend=kk == 3, color=[color_[i] for i in [1, 0]]
        )
        ax[kk].set_ylim([0, 10])
        ax[kk].set_xlim([0, 150])

    infection_rate = 0.15
    vac_init_perc = 0.8
    s_init_perc = 1 - vac_init_perc - inf_init_perc - inf_vac_init_perc
    vac_protect, ef_names = [0.7, 0.8, 0.9], ["Baja", "Moderada", "Alta"]
    results = {}
    for vac_prt in vac_protect:
        results[int(100 * vac_prt)] = sirvid_model(
            days,
            kappa,
            n_population,
            s_init_perc,
            vac_init_perc,
            recovery_rate,
            infection_rate,
            vaccination_rate,
            recovery_rate_vac,
            vac_prt,
            inf_init_perc,
            inf_vac_init_perc,
        )

    infections = pd.concat(
        [aa[0][["infected", "infected_vac"]].sum(axis=1) for ii, aa in results.items()],
        axis=1,
    )
    infections.columns = [
        "{} {}%".format(a, b) for a, b in zip(ef_names, results.keys())
    ]
    infections = 100 * infections / n_population
    ax = infections.iloc[:150].plot(figsize=(6, 4))
    ax.set_ylim([0, 60])
    ax.set_xlim([0, 150])

    deceases = pd.concat(
        [aa[0][["deceased", "deceased_vac"]].sum(axis=1) for ii, aa in results.items()],
        axis=1,
    )
    deceases.columns = ["{} {}%".format(a, b) for a, b in zip(ef_names, results.keys())]
    deceases = 100 * deceases / n_population
    ax = deceases.iloc[:150].plot(figsize=(6, 4), legend=True)
    ax.set_ylim([0, 10])
    ax.set_xlim([0, 150])

    # f, ax = plt.subplots(1, 1)
    # ax.set_xlim([0, max_x])
    # ax.set_ylim([0, 1])
    # plt.show()
    # for i in range(0, max_x + 10, 5):
    #     data_norm.iloc[:i].plot.area(legend=i == 0, ax=ax, color=color_)
    #     plt.pause(0.02)  # pause avec duree en secondes
    #     plt.show()

    """
    SIR
    """
    s_init_perc_sir = 0.99  # S0 (initial susceptible set as % of pop)
    vac_init_perc_sir = 0.0  # V0 (initial vaccinated set as % of pop)
    inf_init_perc_sir = 0.01  # I0 (initial infected set as % of pop)
    inf_vac_init_perc_sir = 0.0  # I0 (initial vaccinated and infected set as % of pop)

    recovery_rate_sir = 1.0  # recovery rate from disease
    vaccination_rate_sir = 0.0  # vaccination rate
    recovery_rate_vac_sir = 0.0  # recovery rate from disease
    vac_protect_sir = 0.0  # decrease of infection rate thanks to the vaccine as %

    assert (
        inf_init_perc_sir + inf_vac_init_perc_sir + s_init_perc_sir + vac_init_perc_sir
        == 1
    )

    data, inc = sirvid_model(
        days,
        kappa,
        n_population,
        s_init_perc_sir,
        vac_init_perc_sir,
        recovery_rate_sir,
        infection_rate,
        vaccination_rate_sir,
        recovery_rate_vac_sir,
        vac_protect_sir,
        inf_init_perc_sir,
        inf_vac_init_perc_sir,
    )

    # charts
    color_sir = [color_[_] for _ in [3, 5, 7]]
    cols = ["susceptibles", "infected", "recovered"]
    data_ = data[cols[::-1]]
    rename_ = {"infected": "infectados", "recovered": "recuperados"}
    data_ = data_.rename(columns=rename_)
    data_norm = data_ / n_population

    f, ax = plt.subplots(1, 2, figsize=(15, 5))
    ax[0].set_xlim([0, max_x])
    ax[0].set_ylim([0, n_population])
    data_.plot(legend=True, ax=ax[0], color=color_sir)

    ax[1].set_xlim([0, max_x])
    ax[1].set_ylim([0, 1])
    data_norm[["recuperados", "infectados", "susceptibles"]].plot.area(
        legend=True, ax=ax[1], color=color_sir
    )
