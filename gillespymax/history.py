from collections import Counter
import h5py
import datetime
import secrets
import polars as pl
from os import PathLike


class ContagionRecords:
    """Records of contagion history.

    Each attribute is a list of values, with one element in each list for each event.
    An empty string is used as an element if the associated field/attribute is not applicable.

    Attributes:
        t : time of event
        enode: node of event that changes state, primary node
        anode: node(s) that are involved but do not change state, auxiliary nodes
        group: group/context through which event occurs
        efrom: primary node's initial state (before event)
        eto: primary node's final state (after event)
        astatus: states of the auxiliary nodes
        txncode: string representation of the states of the nodes involved in the event
    """

    __slots__ = (
        "t",
        "enode",
        "anode",
        "group",
        "efrom",
        "eto",
        "astatus",
        "txncode",
        "states",
    )

    def __init__(self, return_statuses):
        """Initalise an empty contagion record

        Args:
            return_statuses (list): List of states to record a count of at each event
        """
        self.t = []
        self.enode = []
        self.anode = []
        self.group = []
        self.efrom = []
        self.eto = []
        self.astatus = []
        self.txncode = []
        # need to init this way to get distinct lists
        self.states = {state: [] for state in return_statuses}

    def set_initial_condition(self, t0, status):
        """Sets the initial record in the record table

        Args:
            t0 (float) : Initial time
            status (dict) : Initial states (node (str) -> state (str))

        Constructs the first entry in each of the slots. By default, empty strings are used in
        place of Nones.
        """
        self.t.append(t0)
        # append in empty strings for initial condition
        self.enode.append("")
        self.anode.append("")
        self.group.append("")
        self.efrom.append("")
        self.eto.append("")
        self.astatus.append("")
        self.txncode.append("")
        status_counter = Counter(status.values())
        for state, state_record in self.states.items():
            # Counters have a default value of 0
            state_record.append(status_counter[state])

    def add(self, t, enode="", anode="", group="", efrom="", eto="", astatus=""):
        """Add a record of an event to the record table

        Args:
            t (float) : Current time
            enode (str) : Node that is changing state in this event
            anode (str) : Other nodes that affect the event
            group (str) : Group that event occurs in/through (if applicable)
            efrom (str) : State of node prior to change/event
            eto (str) : State of node after change/event
            astatus (str) : status of other nodes that affect the event

        Also generates a txncode that represents the event category
        """
        # compute txn code
        pre_state = ",".join(filter(None, [efrom, astatus]))
        post_state = ",".join(filter(None, [eto, astatus]))
        txncode = f"({pre_state}) -> ({post_state})"

        self.t.append(float(t))
        self.enode.append(str(enode))
        self.anode.append(str(anode))
        self.group.append(str(group))
        self.efrom.append(str(efrom))
        self.eto.append(str(eto))
        self.astatus.append(str(astatus))
        self.txncode.append(txncode)  # already a str
        for state_record in self.states.values():
            state_record.append(state_record[-1])
        self.states.get(efrom, [0])[-1] -= 1
        self.states.get(eto, [0])[-1] += 1

    def write(self, filename: PathLike, sim_id=None, mode="a", **attrs):
        """Output the record table to hdf5 format

        Args:
            filename (str): Path to file to output record table into
            sim_id (Hashable | None, optional): ID to use as the group name in the file. If None, generates a random id. Defaults to None
            mode (str): one of 'a' or 'w'. If 'a', appends; if 'w', overwrites the target file.
            **attrs: Attributes to add to the group
        """

        if sim_id is None:
            sim_id = f"{datetime.datetime.now().timestamp() * 1e6:.0f}_{secrets.token_hex(4)}"

        with h5py.File(filename, mode) as fp:
            grp = fp.create_group(sim_id)
            grp.attrs.update(attrs)
            for attr in (
                "t",
                "enode",
                "anode",
                "group",
                "efrom",
                "eto",
                "astatus",
                "txncode",
            ):
                grp.create_dataset(attr, data=self.__getattribute__(attr))
            for state, record in self.states.items():
                grp.create_dataset(state, data=record)

        return sim_id

    def to_dataframe(self):

        return pl.from_dict(
            {
                **{
                    k: self.__getattribute__(k)
                    for k in (
                        "t",
                        "enode",
                        "anode",
                        "group",
                        "efrom",
                        "eto",
                        "astatus",
                        "txncode",
                    )
                },
                **self.states,
            }
        )


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
