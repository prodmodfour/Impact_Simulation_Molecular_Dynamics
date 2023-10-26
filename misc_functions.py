def flatten_two_d_array(two_d_array):
    flattened_array = list()
    for array in two_d_array:
        for value in array:
            flattened_array.append(value)

    return flattened_array

def remove_duplicates(array):
    new_array = list()

    for value in array:
        if value in new_array:
            continue
        elif value not in new_array:
            new_array.append(value)

    return new_array

def check_if_integer(*args):
    for argument in (args):
        try:
            argument = int(argument)
        except:
            print("Arguments must be integers.")

def print_debug_key(settings, debug_key):
    if settings.debug_keys is True:
        print(debug_key)