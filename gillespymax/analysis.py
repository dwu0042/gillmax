import h5py
import polars as pl
from os import PathLike
from typing import Iterable

def read_records(filename: PathLike):

    dataframes = dict()
    with h5py.File(filename, "r") as fp:
        for group in fp:
            records = {
                k: (list(v.asstr()) if h5py.check_string_dtype(v.dtype) else v[:])
                for k, v in fp[group].items()
            }
            dataframes[group] = pl.from_dict(records)

    return dataframes

def plot_history(dataframe: pl.DataFrame, states_to_plot: Iterable):

    return (
        dataframe
        .select('t', list(states_to_plot))
        .unpivot(
            index='t',
            variable_name='state',
            value_name='count',
        )
        .hvplot.step(
            x='t',
            y='count',
            by='state'
        )
    )