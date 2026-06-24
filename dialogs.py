import tkinter as tk


def ask_user_parameter(parent,
                       title: str,
                       message: str,
                       input_type: str = "int",
                       default=None):
    """
    Reusable modal dialog for requesting a single parameter from the user.

    Args:
        parent: Tkinter parent window.
        title: Dialog window title.
        message: Prompt message displayed to user.
        input_type: Expected type — "int", "float", or "str".
        default: Default value pre-filled in entry field.

    Returns:
        The user's input (cast to input_type), or None if skipped/cancelled.
    """
    result = [None]  # Mutable container for closure

    dialog = tk.Toplevel(parent)
    dialog.title(title)
    dialog.geometry("450x220")
    dialog.transient(parent)
    dialog.grab_set()

    # Message
    tk.Label(dialog, text=message, wraplength=400, justify="left",
             font=("Arial", 10)).pack(pady=(15, 10), padx=15)

    # Entry
    entry_var = tk.StringVar(value=str(default) if default is not None else "")
    entry = tk.Entry(dialog, textvariable=entry_var, width=15, font=("Arial", 11))
    entry.pack(pady=5)
    entry.focus_set()

    # Validation feedback
    error_label = tk.Label(dialog, text="", fg="red", font=("Arial", 9))
    error_label.pack()

    def on_confirm():
        raw = entry_var.get().strip()
        if not raw:
            error_label.config(text="Please enter a value.")
            return
        try:
            if input_type == "int":
                result[0] = int(raw)
            elif input_type == "float":
                result[0] = float(raw)
            else:
                result[0] = raw
            dialog.destroy()
        except ValueError:
            error_label.config(text=f"Invalid input. Expected a {input_type} value.")

    def on_skip():
        dialog.destroy()

    # Buttons
    btn_frame = tk.Frame(dialog)
    btn_frame.pack(pady=15)

    tk.Button(btn_frame, text="Confirm", command=on_confirm,
              padx=10).pack(side="left", padx=10)
    tk.Button(btn_frame, text="Skip", command=on_skip,
              padx=10).pack(side="left", padx=10)

    # Allow Enter key to confirm
    entry.bind("<Return>", lambda e: on_confirm())

    parent.wait_window(dialog)
    return result[0]


def ask_ligand_choice(parent,
                      ligand_1_name: str,
                      ligand_2_name: str,
                      plate_info: str = ""):
    """
    Modal dialog asking the user which of two ligands was used on a specific plate.
    Shown when the protocol defines two ligands but the layout indicates
    only one ligand per plate (so the software cannot infer which one).
    Returns:
        "L1" if Ligand 1 selected, "L2" if Ligand 2 selected, or None if skipped.
    """
    result = [None]

    dialog = tk.Toplevel(parent)
    dialog.title("Select Ligand")
    dialog.geometry("450x250")
    dialog.transient(parent)
    dialog.grab_set()

    plate_str = f" for plate '{plate_info}'" if plate_info else ""
    msg = (f"Protocol defines two ligands but layout is 'one ligand per plate'.\n"
           f"Which ligand was used{plate_str}?")
    tk.Label(dialog, text=msg, wraplength=400, justify="left",
             font=("Arial", 10)).pack(pady=(15, 10), padx=15)

    # Preselect L1
    choice_var = tk.StringVar(value="L1")

    rb_frame = tk.Frame(dialog)
    rb_frame.pack(pady=5, padx=30, anchor="w")

    tk.Radiobutton(rb_frame, text=f"Ligand 1:  {ligand_1_name}",
                   variable=choice_var, value="L1",
                   font=("Arial", 10)).pack(anchor="w", pady=2)
    tk.Radiobutton(rb_frame, text=f"Ligand 2:  {ligand_2_name}",
                   variable=choice_var, value="L2",
                   font=("Arial", 10)).pack(anchor="w", pady=2)

    def on_confirm():
        result[0] = choice_var.get()
        dialog.destroy()

    def on_skip():
        dialog.destroy()

    btn_frame = tk.Frame(dialog)
    btn_frame.pack(pady=15)

    tk.Button(btn_frame, text="Confirm", command=on_confirm,
              padx=10).pack(side="left", padx=10)
    tk.Button(btn_frame, text="Skip", command=on_skip,
              padx=10).pack(side="left", padx=10)

    parent.wait_window(dialog)
    return result[0]


def ask_ligand_layout(parent,
                      ligand_1_name: str,
                      ligand_2_name: str,
                      protocol_name: str = ""):
    """
    Modal dialog asking the user to choose a ligand layout for all plates.
    Shown when the protocol defines two ligands but uses a 'one cell line per plate' cell
    layout without specifying a ligand layout (older protocols).
    Returns:
        "half", "alternating", or "one ligand" if confirmed, or None if skipped.
    """
    result = [None]

    dialog = tk.Toplevel(parent)
    dialog.title("Select Ligand Layout")
    dialog.geometry("500x300")
    dialog.transient(parent)
    dialog.grab_set()

    # Prompt message
    proto_str = f" (Protocol: '{protocol_name}')" if protocol_name else ""
    msg = (f"Protocol defines two ligands but no ligand layout is specified{proto_str}.\n"
           f"How are '{ligand_1_name}' and '{ligand_2_name}' distributed across the plate?\n\n"
           f"This choice will apply to all plates of this protocol.")
    tk.Label(dialog, text=msg, wraplength=450, justify="left",
             font=("Arial", 10)).pack(pady=(15, 10), padx=15)

    # Preselect "half/half"
    choice_var = tk.StringVar(value="half")

    rb_frame = tk.Frame(dialog)
    rb_frame.pack(pady=5, padx=30, anchor="w")

    tk.Radiobutton(rb_frame, text="half/half",
                   variable=choice_var, value="half",
                   font=("Arial", 10)).pack(anchor="w", pady=2)
    tk.Radiobutton(rb_frame, text="alternating",
                   variable=choice_var, value="alternating",
                   font=("Arial", 10)).pack(anchor="w", pady=2)
    tk.Radiobutton(rb_frame, text="one ligand per plate  ( — will be asked to specify per plate)",
                   variable=choice_var, value="one ligand",
                   font=("Arial", 10)).pack(anchor="w", pady=2)

    def on_confirm():
        result[0] = choice_var.get()
        dialog.destroy()

    def on_skip():
        dialog.destroy()

    btn_frame = tk.Frame(dialog)
    btn_frame.pack(pady=15)

    tk.Button(btn_frame, text="Confirm", command=on_confirm,
              padx=10).pack(side="left", padx=10)
    tk.Button(btn_frame, text="Skip", command=on_skip,
              padx=10).pack(side="left", padx=10)

    parent.wait_window(dialog)
    return result[0]

def ask_filename_collision(parent, colliding_files):
    """
    Modal dialog shown for FILENAME_DATA_COLLISION: two or more sources contain a file
    with the SAME name but DIFFERENT raw data. One global choice applies to all colliding
    files — the first occurrence keeps its name, the others are renamed.

    Args:
        parent: Tkinter parent window.
        colliding_files: list of file names that collide across sources.

    Returns:
        "rename" to keep both copies (renaming later occurrences),
        "cancel" or None to abort the merge.
    """
    result = ["cancel"]  # default to cancel if window is closed without a choice

    dialog = tk.Toplevel(parent)
    dialog.title("Filename Collision")
    dialog.geometry("520x320")
    dialog.transient(parent)
    dialog.grab_set()

    files_preview = ", ".join(str(f) for f in colliding_files[:8])
    if len(colliding_files) > 8:
        files_preview += f", … (+{len(colliding_files) - 8} more)"

    msg = ("Some files share the same name across sources but carry DIFFERENT raw data:\n\n"
           f"{files_preview}\n\n"
           "These are not the same measurement. You can keep both copies by renaming the "
           "later occurrences (the first occurrence keeps its original name), or cancel the "
           "merge to fix the inputs yourself.\n\n"
           "This choice applies to ALL colliding files.")
    tk.Label(dialog, text=msg, wraplength=470, justify="left",
             font=("Arial", 10)).pack(pady=(15, 10), padx=15)

    def on_rename():
        result[0] = "rename"
        dialog.destroy()

    def on_cancel():
        result[0] = "cancel"
        dialog.destroy()

    btn_frame = tk.Frame(dialog)
    btn_frame.pack(pady=15)

    tk.Button(btn_frame, text="Rename & keep both", command=on_rename,
              padx=10).pack(side="left", padx=10)
    tk.Button(btn_frame, text="Cancel merge", command=on_cancel,
              padx=10).pack(side="left", padx=10)

    parent.wait_window(dialog)
    return result[0]


def ask_main_plasmids_selection(parent, ordered_tokens, preselected, context=""):
    """
    Multi-select dialog for choosing the global Main_Plasmids when merging
    sources that encode the same plasmids under different Main_Plasmids/Transfection splits.

    Args:
        parent:         Tkinter parent window.
        ordered_tokens: list of (token, count) ordered highest-occurrence first. `count` is
                        shown as a hint (how many distinct conditions contain the token).
        preselected:    list of tokens to check by default (the declared-main intersection).
        context:        optional explanatory line shown above the list.

    Returns:
        The selected list of tokens (original casing) on confirm, or None on skip/close.
    """
    result = [None]

    dialog = tk.Toplevel(parent)
    dialog.title("Select Main Plasmids")
    dialog.geometry("520x420")
    dialog.transient(parent)
    dialog.grab_set()

    msg = ("Select the plasmid(s) that form the shared Main Plasmids for the "
           "merged data. The remaining plasmids become each condition's Transfection.")
    if context:
        msg = context + "\n\n" + msg
    tk.Label(dialog, text=msg, wraplength=480, justify="left",
             font=("Arial", 10)).pack(pady=(15, 10), padx=15, side="top")

    token_vars = []  # list of (token, BooleanVar)

    def on_confirm():
        chosen = [tok for tok, var in token_vars if var.get()]
        if not chosen:
            error_label.config(text="Select at least one plasmid for the backbone.")
            return
        result[0] = chosen
        dialog.destroy()

    def on_skip():
        result[0] = None
        dialog.destroy()

    # Buttons
    btn_frame = tk.Frame(dialog)
    btn_frame.pack(side="bottom", pady=15)
    tk.Button(btn_frame, text="Confirm", command=on_confirm,
              padx=10).pack(side="left", padx=10)
    tk.Button(btn_frame, text="Cancel", command=on_skip,
              padx=10).pack(side="left", padx=10)

    error_label = tk.Label(dialog, text="", fg="red", font=("Arial", 9))
    error_label.pack(side="bottom")

    # Scrollable checkbutton list (one token per row, highest-occurrence first)
    list_frame = tk.Frame(dialog)
    list_frame.pack(side="top", fill="both", expand=True, padx=20, pady=5)
    canvas = tk.Canvas(list_frame, highlightthickness=0)
    scrollbar = tk.Scrollbar(list_frame, orient="vertical", command=canvas.yview)
    inner = tk.Frame(canvas)
    inner.bind("<Configure>",
               lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
    canvas.create_window((0, 0), window=inner, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

    preselected_cf = {str(t).casefold() for t in (preselected or [])}
    for token, count in ordered_tokens:
        var = tk.BooleanVar(value=str(token).casefold() in preselected_cf)
        label = f"{token}   ({count} measurement{'s' if count != 1 else ''})"
        tk.Checkbutton(inner, text=label, variable=var,
                       font=("Arial", 10), anchor="w").pack(anchor="w", fill="x")
        token_vars.append((token, var))

    parent.wait_window(dialog)
    return result[0]


def warn_and_abort(parent, title, message):
    """
    Simple warning dialog with no proceed option. Generic sink for every no-override
    forbidden merge code (TIME_VECTOR_DIVERGENCE, LABELING_MISMATCH, RAW_CHANNELS_REQUIRED,
    and any future no-override code). The caller passes the issue's own code as `title`
    and its own `message`; this dialog does not author per-code text.

    Returns:
        None (informational only).
    """
    dialog = tk.Toplevel(parent)
    dialog.title(str(title) if title else "Merge Blocked")
    dialog.geometry("520x280")
    dialog.transient(parent)
    dialog.grab_set()

    tk.Label(dialog, text=str(title), fg="#b00000",
             font=("Arial", 11, "bold")).pack(pady=(15, 5), padx=15)
    tk.Label(dialog, text=str(message), wraplength=470, justify="left",
             font=("Arial", 10)).pack(pady=(0, 10), padx=15)

    tk.Button(dialog, text="OK", command=dialog.destroy,
              padx=20).pack(pady=15)

    parent.wait_window(dialog)
    return None