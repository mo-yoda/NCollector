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