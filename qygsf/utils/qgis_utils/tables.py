
import numbers


def field_values(
    layer,
    curr_field_ndx: numbers.Integral
):

    values = []

    if layer.selectedFeatureCount() > 0:
        features = layer.selectedFeatures()
    else:
        features = layer.getFeatures()

    for feature in features:
        values.append(feature.attributes()[curr_field_ndx])

    return values


def vect_attrs(
    layer,
    field_list
):

    if layer.selectedFeatureCount() > 0:
        features = layer.selectedFeatures()
    else:
        features = layer.getFeatures()

    provider = layer.dataProvider()
    field_indices = [provider.fieldNameIndex(field_name) for field_name in field_list]

    # retrieve (selected) attributes features
    data_list = []
    for feature in features:
        attrs = feature.fields().toList()
        data_list.append([feature.attribute(attrs[field_ndx].name()) for field_ndx in field_indices])

    return data_list

