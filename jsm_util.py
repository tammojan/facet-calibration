from __future__ import print_function

def get_config(parameter, config, source=None, default=None):
    """
    Gets the value of a parameter from the config variable.
    It looks in the following order:
      * a specific value that depends on the source
      * a default value that does not depends on the source
      * the default parameter entered in 'default'
    """
    if source is None:
        try:
            return config[parameter]["default"]
        except (KeyError, TypeError):
            if default is not None:
                return config.get(parameter, default)
            else:
                return config[parameter]
    else:
        if parameter in config.keys():
            if "default" in config[parameter].keys():
                return config[parameter].get(source, config[parameter]["default"])
            else:
                if default is not None:
                    return config[parameter].get(source, default)
                else:
                    return config[parameter][source]
        else:
            if default is not None:
                return config.get(parameter, default)
            else:
                return config[parameter]
            
