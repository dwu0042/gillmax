import h5py
import polars as pl
import seaborn as sns
from os import PathLike
from typing import Iterable

def read_records(filename: PathLike):
    """Returns a Mapping of simulation ids to parsed simulation dataframes"""
    dataframes = dict()
    with h5py.File(filename, "r") as fp:
        for group in fp:
            records = {
                k: (list(v.asstr()) if h5py.check_string_dtype(v.dtype) else v[:])
                for k, v in fp[group].items()
            }
            dataframes[group] = pl.from_dict(records)

    return dataframes

def _df_longer(dataframe: pl.DataFrame, states_to_plot: Iterable):
    return (
        dataframe
        .select('t', *list(states_to_plot))
        .unpivot(
            index='t',
            variable_name='state',
            value_name='count',
        )
    )

def _plot_hist(df_long: pl.DataFrame, ax=None, **kwargs):
    return sns.lineplot(
        df_long,
        x='t',
        y='count',
        hue='state',
        drawstyle='steps',
        ax=ax,
        **kwargs
    )

def plot_history(dataframe: pl.DataFrame, states_to_plot: Iterable, ax=None, **kwargs):
    """Plots the state history for the given states to plot"""

    return _plot_hist(_df_longer(dataframe, states_to_plot), ax=ax, **kwargs)

def compute_delay_distribution(dataframe: pl.DataFrame, init_code, completion_code):
    """Computes the time delay between two types of events, where continuity is provided by being the same individual (enode)"""

    init_enodes = dataframe.filter(pl.col('txncode').eq(init_code)).select('enode').to_series()
    complete_enodes = dataframe.filter(pl.col('txncode').eq(completion_code)).select('enode').to_series()

    return (
        dataframe
        .filter(
            pl.col('enode').is_in(init_enodes),
            pl.col('enode').is_in(complete_enodes),
        )
        .select(
            't', 'enode', 'txncode',
        )
        .pivot(index='enode', on='txncode')
        .select(
            'enode',
            (
                pl.col(completion_code) - pl.col(init_code)
            ).alias('delay')
        )
    )