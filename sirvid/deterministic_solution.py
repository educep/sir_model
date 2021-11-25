"""
Created by ecepeda at 19/11/2021
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

color_ = sns.color_palette("Paired")


def sirvid_model(days, kappa, n_population, s_init_perc, vac_init_perc, recovery_rate,
                 infection_rate, vaccination_rate, recovery_rate_vac, vac_protect):
    """
        defining parameters
    """
    # days = 15  # average duration of disease
    # kappa = 2.5  # avg nb of interactions by day by individual
    # n_population = 1e6  # population size
    # s_init_perc = .985  # S0 (initial susceptible set as % of pop)
    # vac_init_perc = .005  # V0 (initial vaccinated set as % of pop)
    inf_vac_init_perc = vac_init_perc * 0.001
    inf_init_perc = 1 - s_init_perc - vac_init_perc - inf_vac_init_perc  # I0 (initial infected set as % of pop)
    assert inf_init_perc + inf_vac_init_perc + s_init_perc + vac_init_perc == 1

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
    # vac_protect = .95  # decrease of infection rate thanks to vaccin as %
    infection_rate_vac = infection_rate * (1 - vac_protect)  # infection proba between vaccinated and infected
    beta_vac = kappa * infection_rate_vac  # transmission rate for vaccinated people
    gamma_vac = recovery_rate_vac / days  # recovery rate by day
    mu_vac = mortality_vac / days  # decease rate by day
    nbdays = 150  # total number of days

    # declaring objects
    time_frame = range(0, nbdays)
    data = pd.DataFrame(0, index=time_frame, columns=['susceptibles', 'infected', 'recovered', 'deceased', 'vaccinated',
                                                      'infected_vac', 'recovered_vac', 'deceased_vac'])
    increments = pd.DataFrame(0, index=time_frame, columns=['vaccinations', 'infections', 'recoveries', 'deceases',
                                                            'infections_vac', 'recoveries_vac', 'deceases_vac'])

    # initializing instances, the other population compartments are set to 0 (always)
    # problem, if vaccinated people are significant, say > 10%, we should consider a % > 0 of infected_vac
    data.loc[0, 'susceptibles'] = n_population * s_init_perc
    data.loc[0, 'vaccinated'] = n_population * vac_init_perc
    data.loc[0, 'infected'] = n_population * inf_init_perc
    data.loc[0, 'infected_vac'] = n_population * inf_vac_init_perc

    for day_ in range(1, nbdays):
        # day_ = 1
        prev_day = day_ - 1
        # computing increments for non vaccinated people
        suscept_prev_day = data.loc[prev_day, 'susceptibles']
        infected_prev_day = data.loc[prev_day, 'infected']
        total_infected_prev_day = infected_prev_day + data.loc[prev_day, 'infected_vac']
        alive_pop_prev_day = data.loc[prev_day, ['susceptibles', 'infected', 'recovered',
                                                 'vaccinated', 'infected_vac', 'recovered_vac']].sum()

        increments.loc[day_, 'vaccinations'] = vaccination_rate * suscept_prev_day
        increments.loc[day_, 'infections'] = beta * suscept_prev_day * total_infected_prev_day / alive_pop_prev_day
        increments.loc[day_, 'recoveries'] = gamma * infected_prev_day
        increments.loc[day_, 'deceases'] = mu * infected_prev_day

        # computing increments for vaccinated people
        vac_prev_day = data.loc[prev_day, 'vaccinated']
        infected_vac_prev_day = data.loc[prev_day, 'infected_vac']
        increments.loc[day_, 'infections_vac'] = beta_vac * vac_prev_day * total_infected_prev_day / alive_pop_prev_day
        increments.loc[day_, 'recoveries_vac'] = gamma_vac * infected_vac_prev_day
        increments.loc[day_, 'deceases_vac'] = mu_vac * infected_vac_prev_day

        # aggregating increments for non vaccinated people
        data.loc[day_, 'susceptibles'] = suscept_prev_day - increments.loc[day_, ['vaccinations', 'infections']].sum()
        data.loc[day_, 'infected'] = infected_prev_day + increments.loc[day_, 'infections'] - \
                                     increments.loc[day_, ['recoveries', 'deceases']].sum()
        data.loc[day_, 'recovered'] = data.loc[prev_day, 'recovered'] + increments.loc[day_, 'recoveries']
        data.loc[day_, 'deceased'] = data.loc[prev_day, 'deceased'] + increments.loc[day_, 'deceases']

        # aggregating increments for vaccinated people
        data.loc[day_, 'vaccinated'] = data.loc[prev_day, 'vaccinated'] + increments.loc[day_, 'vaccinations'] - \
                                       increments.loc[day_, 'infections_vac']
        data.loc[day_, 'infected_vac'] = data.loc[prev_day, 'infected_vac'] + increments.loc[day_, 'infections_vac'] \
                                         - increments.loc[day_, ['recoveries_vac', 'deceases_vac']].sum()
        data.loc[day_, 'recovered_vac'] = data.loc[prev_day, 'recovered_vac'] + increments.loc[day_, 'recoveries_vac']
        data.loc[day_, 'deceased_vac'] = data.loc[prev_day, 'deceased_vac'] + increments.loc[day_, 'deceases_vac']

    return data


if __name__ == '__main__':
    """
    defining parameters
    """
    days = 15  # average duration of disease
    kappa = 2.5  # avg nb of interactions by day by individual
    n_population = 1e6  # population size
    s_init_perc = .985  # S0 (initial susceptible set as % of pop)
    vac_init_perc = .005  # V0 (initial vaccinated set as % of pop)
    # inf_vac_init_perc = vac_init_perc * 0.001
    # inf_init_perc = 1 - s_init_perc - vac_init_perc - inf_vac_init_perc  # I0 (initial infected set as % of pop)
    # assert inf_init_perc + inf_vac_init_perc + s_init_perc + vac_init_perc == 1
    #
    # parameters for non vaccinated people
    recovery_rate = .90  # recovery rate from disease
    mortality = 1 - recovery_rate  # rate decease from disease
    infection_rate = 0.14  # probability of infection for one interaction between susceptible and infected
    # beta = kappa * infection_rate  # transmission rate
    # gamma = recovery_rate / days  # recovery rate by day
    # mu = mortality / days  # decease rate by day
    #
    # parameters for vaccinated people
    vaccination_rate = 0.005  # vaccination rate
    recovery_rate_vac = .99  # recovery rate from disease
    # mortality_vac = 1 - recovery_rate_vac  # rate decease from disease
    vac_protect = .95  # decrease of infection rate thanks to vaccin as %
    # infection_rate_vac = infection_rate * (1 - vac_protect)  # infection proba between vaccinated and infected
    # beta_vac = kappa * infection_rate_vac  # transmission rate for vaccinated people
    # gamma_vac = recovery_rate_vac / days  # recovery rate by day
    # mu_vac = mortality_vac / days  # decease rate by day
    # nbdays = 150  # total number of days

    data = sirvid_model(days, kappa, n_population, s_init_perc, vac_init_perc, recovery_rate,
                        infection_rate, vaccination_rate, recovery_rate_vac, vac_protect)

    # charts
    cols = ['susceptibles', 'vaccinated', 'infected', 'infected_vac', 'recovered', 'recovered_vac',
            'deceased', 'deceased_vac']
    data_norm = data[cols[::-1]] / n_population
    # data.iloc[:100].plot.area()
    max_x = 80
    f, ax = plt.subplots(1, 1)
    ax.set_xlim([0, max_x])
    ax.set_ylim([0, 1])
    plt.show()
    for i in range(0, max_x + 10, 5):
        data_norm.iloc[:i].plot.area(legend=i == 0, ax=ax, color=color_)
        plt.pause(0.02)  # pause avec duree en secondes


