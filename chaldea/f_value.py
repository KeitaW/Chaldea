import numpy as np
def take_correspondence(labels_true, labels_pred):
    """
    -1:(noise)も扱えるように改良したバージョン
    """
    clust_labels = np.unique(labels_true)
    corresponding_labels_pred = np.ones_like(labels_pred) * (-1)
    correspondence_table = dict()

    for label in clust_labels[clust_labels != -1]:
        predicted_labels, labels_counts = np.unique(labels_pred[labels_true == label], return_counts=True)
        correspondence_table[label] = predicted_labels[np.argmin(np.abs(labels_counts-24))]
        #correspondence_table[label] = predicted_labels[np.argmax(labels_counts)]
        corresponding_labels_pred[labels_pred == correspondence_table[label]] = label
    print(correspondence_table)
    for label in np.unique(labels_pred):
        if label in correspondence_table.values():
            continue
        else:
            corresponding_labels_pred[labels_pred == label] = -1            
    return(corresponding_labels_pred)
def purity(labels_true, labels_pred):
    labels_pred_ = take_correspondence(labels_true, labels_pred)
    return(np.sum([np.sum(labels_pred_[labels_true == label] == label) for label in np.unique(labels_true)])/len(labels_true))
def inv_purity(labels_true, labels_pred):
    labels_pred_ = take_correspondence(labels_true, labels_pred)
    ipurity = 0
    for label in np.unique(labels_true):
        ipurity += np.sum(labels_pred_ == label)* np.sum(np.logical_and(labels_true == label, labels_pred_ == label)) / np.sum(labels_true == label) 
    ipurity /= len(labels_true)
    return(ipurity)
def f_value(labels_true, labels_pred):
    print("labels")
    print(labels_true[:100])
    print(labels_pred[:100])
    purity_ = purity(labels_true, labels_pred)
    print(purity_)
    ipurity_ = inv_purity(labels_true, labels_pred)
    print(ipurity_)
    return(2*purity_*ipurity_/(purity_+ipurity_))
