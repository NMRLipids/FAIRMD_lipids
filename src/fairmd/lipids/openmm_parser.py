from xml.etree import ElementTree as etree


def openmm_input(input, type):
    if type == "inp":
        return read_inp_input(input)

    if type == "xml":
        return read_xml_input(input)


def read_inp_input(input):
    try:
        filename = open(input)
        print(filename.name)

    except FileNotFoundError:
        print("File does not exist")

    with open(filename.name) as f:
        tmp_type = []
        tmp_value = []

        for line in f:
            split = line.split()

            if len(split) > 2:
                tmp_type.append(split[0])
                tmp_value.append(split[2])

    return [tmp_type, tmp_value]


def read_xml_input(input):
    try:
        filename = open(input)

    except FileNotFoundError:
        print("File does not exist")

    tree = etree.parse(filename.name)
    temp = [value for key, value in tree.getroot()[1].attrib.items() if "temperature" in key.lower()]

    return (["temp"], [temp])


class OpenMMParser:
    def __init__(self, filename, type):
        [properties, values] = openmm_input(filename, type)
        print(properties)
        self.temperature = float(values[properties.index("temp")])

    #       self.pressure = values[properties.index('p_ref')]
    #       self.ptype = values[properties.index('p_type')]
    #       self.dt = values[properties.index('dt')]

    def set_tempurature(self, temp):
        self.temperature = temp

    def set_pressure(self, pressure):
        self.pressure = pressure

    def set_ptype(self, ptype):
        self.ptype = ptype

    def set_dt(self, dt):
        self.dt = dt
