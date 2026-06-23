"""Tkinter-free helpers for main-plasmids grouping and selection."""

NO_COMMON_PLASMIDS = "No Common Plasmids Detected"


def group_by_main_plasmids(experiment):
    """Group folders by ``tuple(folder.protocol.main_plasmids)``.

    Folders without a protocol, or whose protocol has no (falsy) ``main_plasmids``,
    are skipped. Insertion order is preserved (the first set encountered stays first),
    matching the historical "default to the first set" behaviour.

    Args:
        experiment: iterable of MeasurementFolder-like objects. Only
            ``folder.protocol.main_plasmids`` is read.

    Returns:
        dict mapping ``tuple(main_plasmids)`` -> list of folders for that set.
    """
    groups: dict[tuple, list] = {}
    for folder in experiment:
        protocol = getattr(folder, "protocol", None)
        if protocol is None:
            continue
        main_plasmids = getattr(protocol, "main_plasmids", None)
        if not main_plasmids:
            continue
        key = tuple(main_plasmids)
        groups.setdefault(key, []).append(folder)
    return groups


def resolve_selection(experiment, chooser):
    """Group ``experiment`` and resolve which main-plasmids set to keep.

    This is the core shared by the Import tab and the Merge tab. The
    interactive part (picking among multiple sets) is delegated to ``chooser``,
    which the caller injects (a Tk dialog in the GUI).

    Args:
        experiment: list of MeasurementFolder-like objects.
        chooser: callable ``chooser(groups: dict) -> key | None`` invoked **only**
            when more than one set is present. It receives the ordered
            ``{tuple: [folders]}`` mapping and returns the chosen tuple key, or
            ``None`` to signal cancel/abort.

    Returns:
        ``(filtered_experiment, selected_str)`` on success:
          * no set declared      -> ``(experiment, "No Common Plasmids Detected")``
          * exactly one set      -> ``(experiment, "<a + b>")`` (chooser NOT called)
          * multiple, a chosen   -> ``(folders_for_chosen_set, "<a + b>")``
        ``(None, None)`` if the user cancels the multi-set prompt.

    The passed ``experiment`` is never mutated. For the single/no-set cases the
    same list object is returned unchanged so callers can preserve every folder
    (including protocol-less ones) exactly as before.
    """
    groups = group_by_main_plasmids(experiment)

    if not groups:
        return experiment, NO_COMMON_PLASMIDS

    if len(groups) == 1:
        only_key = next(iter(groups))
        return experiment, " + ".join(only_key)

    # More than one set: ask the caller to pick one.
    chosen_key = chooser(groups)
    if chosen_key is None:
        return None, None

    return groups[chosen_key], " + ".join(chosen_key)